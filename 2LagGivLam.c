/*
 * =====================================================================================
 *
 *       Filename:  2LagGivLam.c
 *
 *    Description:  When given lambda, this program iterates the value function.
 *
 *        Version:  1.0
 *        Created:  2014年12月15日 13时38分31秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LU Tianshu (), tssslu@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>

#define MAX_PERIOD		100
#define MAX_X			200
#define MAX_D_LENGTH		30
#define N			2
#define INF			99999999
#define MAX(x,y)		((x>y)?x:y)
#define MIN(x,y)		((x<y)?x:y)

//Value Function
double V[MAX_PERIOD][MAX_X*2][2];

//Incremental Holding Cost
double h[N];

//Penalty Cost
double p;

//Demand Distribution and Expected Demand
int D_len;
int D[MAX_D_LENGTH];
double P[MAX_D_LENGTH];

//Capacity Constraint
int K[N];

//Discount Factor
double beta;

//Lower & Upper Bound of X, Number of Periods
int LB, UB, period;

//Lambda
double lambda;

double get_value(int X, int prd, int ech)
{
	return V[prd][X+MAX_X][ech];
}

void set_value(int X, int prd, int ech, int val)
{
	V[prd][X+MAX_X][ech] = val;
}

double l(int Y)
{
	int i;
	double tmpEV = 0;
	for (i = 0; i < D_len; i++) {
		if (Y >= D[i]) {
			tmpEV += (h[0]+h[1]) * (Y-D[i]) * P[i];
		}
		else {
			tmpEV += p * (D[i]-Y) * P[i];
		}
	}
	return tmpEV;
}

double J1(int Y1, int prd)
{
	int i;
	double res, tmpEV = 0;
	res = l(Y1) - h[1] * Y1 + lambda * Y1;
	for (i = 0; i < D_len; i++) {
		tmpEV += P[i] * get_value(Y1-D[i], prd-1, 0);
	}
	return res + beta * tmpEV;
}

double J2(int Y2, int prd)
{
	int i;
	double tmpEV = 0;
	for (i = 0; i < D_len; i++) {
		tmpEV += P[i] * get_value(Y2-D[i], prd-1, 1);
	}
	return h[1] * Y2 + beta * tmpEV;
}

void DP1(int prd)
{
	int X1, Y1;
	double tmpmin = INF, tmpJ;
	for (X1 = LB; X1 <= UB; X1 ++) {
		for (Y1 = X1; 1; Y1 ++) {
			tmpJ = J1(Y1, prd);
			if (tmpJ < tmpmin) {
				tmpmin = tmpJ;
			}
			else {
				break;
			}
		}
		set_value(X1, prd, 0, tmpmin - lambda * (X1 + K[0]));
	}
}

void DP2(int prd)
{
	int X2, Y1;
	double tmpmin1, tmpmin2, tmpJ;
	tmpmin1 = tmpmin2 = INF;
	for (X2 = LB; X2 <= UB; X2 ++) {
		for (Y1 = X2; 1; Y1 ++) {
			tmpJ = J1(Y1, prd);
			if (tmpJ < tmpmin1) {
				tmpmin1 = tmpJ;
			}
			else {
				break;
			}
		}
		for (Y1 = X2; Y1 <= X2 + K[1]; Y1++) {
			tmpJ = J2(Y1, prd);
			if (tmpJ < tmpmin2) {
				tmpmin2 = tmpJ;
			}
			else {
				break;
			}
		}
		set_value(X2, prd, 1, J1(X2, prd)-tmpmin1+tmpmin2);
	}
}

void DP()
{
	int prd;
	for (prd = 1; prd <= period; prd++) {
		DP1(prd);
		DP2(prd);
	}
}

double J(int Y[], int prd)
{
	int i;
	double tmpEV = 0;
	for (i = 0; i < D_len; i++) {
		tmpEV += P[i] * (get_value(Y[0]-D[i], prd-1, 0)+
				get_value(Y[1]-D[i], prd-1, 1));
	}
	return l(Y[0]) + h[1]*(Y[1]-Y[0]) + beta * tmpEV;
}

void get_policy(int X[], int prd, int plc[])
{
	int Y[2], YUB[2];
	double tmpJ, tmpmin = INF;
	YUB[0] = MIN(X[0], X[0] + K[0]);
	YUB[1] = X[1] + K[1];
	for (Y[0] = X[0]; Y[0] <= YUB[0]; Y[0] ++) {
		for (Y[1] = X[1]; Y[1] <= YUB[1]; Y[1] ++) {
			tmpJ = J(Y, prd);
			if (tmpJ < tmpmin) {
				tmpmin = tmpJ;
				plc[0] = Y[0];
				plc[1] = Y[1];
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
	for (i = 0; i < MAX_PERIOD; i++) {
		for (j = 0; j < MAX_X *2; j++) {
			for (k = 0; k < 2; k++) {
				V[i][j][k] = 0;
			}
		}
	}
}

int main(int argc, const char *argv[])
{
	int X[2], Plc[2];
	scanf("%lf", &lambda);
	init();
	DP();
	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		get_policy(X, period, Plc);
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2lf\n", X[0], X[1]-X[0], X[0], X[1], 
				Plc[0], Plc[1], Plc[0]-X[0], Plc[1]-X[1],
				get_value(X[0], period, 0)+get_value(X[1],period,1));
	}
	return 0;
}
