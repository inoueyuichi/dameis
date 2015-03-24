/*
 * =====================================================================================
 *
 *       Filename:  2Stage.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2015年01月21日 14时23分43秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LU Tianshu (), tssslu@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdio.h>

#define MAX_PERIOD		100
#define MAX_X			200
#define MAX_D_LENGTH		30
#define N			2
#define INF			99999999
#define MAX(x,y)		((x>y)?x:y)
#define MIN(x,y)		((x<y)?x:y)

//Value Function
double V[MAX_PERIOD][MAX_X*2][2];
double SV[MAX_PERIOD][MAX_X*2][MAX_X*2];

//Policy
int Plc[MAX_PERIOD][MAX_X*2][2];
int SPlc[MAX_PERIOD][MAX_X*2][MAX_X*2];

//Ending Echelon Inventory
int Ending[2];

//Incremental Holding Cost
double h[N];

//Penalty Cost
double p;

//Demand Distribution and Expected Demand
int D_len;
int D[MAX_D_LENGTH];
double P[MAX_D_LENGTH];
double ED = 0;

//Capacity Constraint
int K[N];

//Discount Factor
double beta;

//Lower & Upper Bound of X, Number of Periods
int LB, UB, period;

//Initial Installation & Echelon Inventory
int x[2], X[2];

//Langrange Factor
double lambda;

double get_value(int prd, int ech, int X)
{
	return V[prd][MAX_X+X][ech];
}

void set_value(int prd, int ech, int X, int val)
{
	V[prd][MAX_X+X][ech] = val;
}

double get_simu_value(int prd, int X[])
{
	return SV[prd][MAX_X+X[0]][MAX_X+X[1]];
}

void set_simu_value(int prd, int X[], int val)
{
	SV[prd][MAX_X+X[0]][MAX_X+X[1]] = val;
}

void set_policy(int prd, int ech, int X, int plc)
{
	Plc[prd][MAX_X+X][ech] = plc;
}

double l(int Y)
{
	int i;
	double res;
	for (i = 0; i < D_len; i++) {
		if (Y >= D[i]) {
			res += P[i] * (h[1]+h[2]) * (Y-D[i]);
		}
		else {
			res += P[i] * p * (D[i] - Y);
		}
	}
	return res;
}

double J1(int Y1, int prd)
{
	int i;
	double res, tmp = 0;
	res = l(Y1) - h[1] * Y1;
	for (i = 0; i < D_len; i++) {
		tmp += P[i] * get_value(prd-1, 0, Y1-D[i]);
	}
	return res + beta * tmp;
}

double J2(int Y2, int prd)
{
	int i;
	double tmp = 0;
	for (i = 0; i < D_len; i++) {
		tmp += P[i] * get_value(prd-1, 1, Y2-D[i]);
	}
	return beta * tmp + h[1] * Y2;
}

double E1_S1(int X, int prd, int * resY)
{
	int Y1;
	double tmp, tmpmin = INF;
	for (Y1 = X; Y1 <= UB; Y1 ++) {
		tmp = J1(Y1, prd);
		if (tmp< tmpmin) {
			tmpmin = tmp;
			*resY = Y1;
		}
		else {
			break;
		}
	}
	return tmpmin;
}

double E1_S2(int X, int prd, int * resY)
{
	int Y1;
	double tmp, tmpmin = INF;
	for (Y1 = X; Y1 <= UB; Y1 ++) {
		tmp = J1(Y1, prd) + lambda * Y1;
		if (tmp< tmpmin) {
			tmpmin = tmp;
			*resY = Y1;
		}
		else {
			break;
		}
	}
	return tmpmin - lambda * (X + K[0]);
}

double E2(int X2, int prd, int * resY)
{
	int Y2;
	double tmp, tmpmin = INF, res;
	for (Y2 = X2; Y2 <= X2+ K[1]; Y2 ++) {
		tmp = J2(Y2, prd);
		if (tmp < tmpmin) {
			tmpmin = tmp;
			*resY = Y2;
		}
		else {
			break;
		}
	}
	res = tmpmin;
	tmpmin = INF;
	for (Y2 = X2; Y2 <= UB; Y2 ++) {
		tmp = J2(Y2, prd);
		if (tmp < tmpmin) {
			tmpmin = tmp;
		}
		else {
			break;
		}
	}
	return res - tmpmin + J1(X2, prd);
}

void get_policy(int prd, int X[], int Y[])
{
	int tmpY[2];
	E1_S1(X[0], prd, &tmpY[0]);
	E2(X[1], prd, &Y[1]);
	if (Y[0] > X[0] + K[0]) {
		E1_S2(X[0], prd, &Y[0]);
	}
	else {
		Y[0] = tmpY[0];
	}
}

void DP()
{
	int prd, X[2], Y[2];
	double tmpV;
	for (prd = 1; prd <= period; prd++) {
		for (X[0] = LB; X[0] <= UB; X[0] ++) {
			tmpV = E1_S1(X[0], prd, &Y[0]);
			if (Y[0] > X[0] + K[0]) {
				tmpV = E1_S2(X[0], prd, &Y[0]);
			}
			set_value(prd, 0, X[0], tmpV);
			set_policy(prd, 0, X[0], Y[0]);
		}
		for (X[1] = LB; X[1] <= UB; X[1] ++) {
			tmpV = E2(X[1], prd, &Y[1]);
			set_value(prd, 1, X[1], tmpV);
			set_policy(prd, 1, X[1], Y[1]);
		}
	}
}

void simu()
{
	int prd, X[2], Y[2], i, tmpY[2];
	double res = 0, tmpEV;
	for (prd = 1; prd <= period; prd ++) {
		for (X[0] = LB; X[0] <= UB; X[0] ++) {
			for (X[1] = LB; X[1] <= UB; X[1] ++) {
				get_policy(prd, X, Y);
				tmpEV = 0;
				for (i = 0; i < D_len; i++) {
					tmpY[0] = Y[0] - D[i];
					tmpY[1] = Y[1] - D[i];
					tmpEV += P[i] * get_simu_value(prd-1, tmpY);
				}
				set_simu_value(prd, X, tmpEV+l(Y[0])+h[1]*(Y[1]-Y[0]));
			}
		}
	}
}

void init()
{
	//TODO: This function initialize all variables including the value
	//	function array.
	//	This part can be changed in order to read data from a file.

	int i, j, k;
	FILE * fp = fopen("2Echelon.dat", "r");
	fscanf(fp, "%lf%d%lf%lf%d%d%lf%d%d%d", &beta,
			&period, &h[0], &h[1], &K[0], &K[1],
			&p, &UB, &LB, &D_len);
	for (i = 0; i < D_len; i++) {
		fscanf(fp, "%d", &D[i]);
	}
	for (i = 0; i < D_len; i++) {
		fscanf(fp, "%lf", &P[i]);
	}
	for (i = 0; i < D_len; i++) {
		ED += D[i] * P[i];
	}
}

int main(int argc, const char *argv[])
{
	int X[2], Y[2];
	init();
	scanf("%lf", &lambda);
	DP();
	simu();
	while (scanf("%d %d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		get_policy(period, X, Y);
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2lf\n",
			X[0], X[1]-X[0], X[0], X[1],
			Y[0], Y[1], Y[0] - X[0], Y[1] - X[1],
			get_simu_value(period, X)
		      );
	}
	return 0;
}
