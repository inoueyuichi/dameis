/*
 * =====================================================================================
 *
 *       Filename:  2Penalty.c
 *
 *    Description:  This algorithm uses penalty functions to approximate the multi-
 		    echelon value function.
 *
 *        Version:  1.0
 *        Created:  2015年03月23日 15时33分12秒
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

#define MAX_PERIOD		35
#define MAX_X			90
#define MAX_D_LENGTH		30
#define INF			99999999
#define MAX(x,y)		((x>y)?x:y)
#define MIN(x,y)		((x<y)?x:y)

//Value Function
double V[MAX_PERIOD][MAX_X*2][3];
double EvV[MAX_PERIOD][MAX_X*2][MAX_X*2][MAX_X*2];

//Policy
int Base[MAX_PERIOD][MAX_X*2][3];
int Plc[MAX_X*2][MAX_X*2][MAX_X*2][3];

//Ending Echelon Inventory
int Ending[3];

//Incremental Holding Cost
double h[3];

//Penalty Cost
double p;

//Demand Distribution and Expected Demand
int D_len;
int D[MAX_D_LENGTH];
double P[MAX_D_LENGTH];

//Capacity Constraint
int K[3];

//Discount Factor
double beta;

//Lower & Upper Bound of X, Number of Periods
int N, LB, UB, period;

void set_EvV(int prd, int X[], double val)
{
	EvV[prd][MAX_X+X[0]][MAX_X+X[1]][MAX_X+X[2]] = val;
}

double get_EvV(int prd, int X[])
{
	return EvV[prd][MAX_X+X[0]][MAX_X+X[1]][MAX_X+X[2]];
}

void set_value(int prd, int X, int cnt, double val)
{
	V[prd][X+MAX_X][cnt] = val;
}

double get_value(int prd, int X, int cnt)
{
	return V[prd][X+MAX_X][cnt];
}

void set_base(int prd, int X, int cnt, int val)
{
	Base[prd][X+MAX_X][cnt] = val;
}

int get_base(int prd, int X, int cnt)
{
	return Base[prd][X+MAX_X][cnt];
}

void set_policy(int X[], int Y[])
{
	Plc[MAX_X+X[0]][MAX_X+X[1]][MAX_X+X[2]][0] = Y[0];
	Plc[MAX_X+X[0]][MAX_X+X[1]][MAX_X+X[2]][1] = Y[1];
	Plc[MAX_X+X[0]][MAX_X+X[1]][MAX_X+X[2]][2] = Y[2];
}

int get_policy(int X[], int cnt)
{
	return Plc[MAX_X+X[0]][MAX_X+X[1]][MAX_X+X[2]][cnt];
}

double l(int Y1)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		if (Y1 > D[i]) {
			res += P[i] * (h[0]+h[1]+h[2]) * (Y1-D[i]);
		}
		else {
			res += P[i] * p * (D[i]-Y1);
		}
	}
	return res;
}

double L1(int Y1, int prd)
{
	int i;
	double EV = 0;
	for (i = 0; i < D_len; i++) {
		EV += P[i] * get_value(prd-1, Y1-D[i], 0);
	}
	return l(Y1) - (h[1]+h[2])*Y1 + beta*EV;
}

double L2(int Y2, int prd)
{
	int i;
	double EV = 0;
	for (i = 0; i < D_len; i++) {
		EV += P[i] * get_value(prd-1, Y2-D[i], 1);
	}
	return h[1]*Y2 + beta*EV;
}

double L3(int Y3, int prd)
{
	int i;
	double EV = 0;
	for (i = 0; i < D_len; i++) {
		EV += P[i] * get_value(prd-1, Y3-D[i], 2);
	}
	return h[2]*Y3 + beta*EV;
}

double Ev_L(int Y[], int prd)
{
	int i, X[3];
	double res = l(Y[0])+h[1]*(Y[1]-Y[0])+h[2]*(Y[2]-Y[0]),
	       EV = 0;
	for (i = 0; i < D_len; i++) {
		X[0] = Y[0] - D[i];
		X[1] = Y[1] - D[i];
		X[2] = Y[2] - D[i];
		EV += P[i] * get_EvV(prd-1, X);
	}
	return res + beta * EV;
}

double PES2(int X2, int prd)
{
	int Ystar, Y1;
	double tmp, tmpmin = INF;
	for (Y1 = LB; Y1 <= UB; Y1++) {
		tmp = L1(Y1, prd);
		if (tmp < tmpmin) {
			tmpmin = tmp;
			Ystar = Y1;
		}
	}
	if (Ystar <= X2) {
		return 0;
	}
	else {
		return L1(X2, prd) - tmpmin;
	}
}

double PES3(int X3, int prd)
{
	int Ystar, Y2;
	double tmp, tmpmin = INF;
	for (Y2 = LB; Y2 < UB; Y2++) {
		tmp = L2(Y2, prd);
		if (tmp < tmpmin) {
			tmpmin = tmp;
			Ystar = Y2;
		}
	}
	if (Ystar <= X3) {
		return 0;
	}
	else {
		return L2(X3, prd) - tmpmin;
	}
}

void DP()
{
	int prd, X, Y, Ystar;
	double tmpmin, tmp;
	for (prd = 1; prd <= period; prd++) {
		if (prd == period) {
			printf(" ");
		}
		for (X = LB; X <= UB; X++) {
			tmpmin = INF;
			for (Y = X; Y <= X+K[0]; Y++) {
				tmp = L1(Y, prd);
				if (tmp < tmpmin) {
					tmpmin = tmp;
					Ystar = Y;
				}
			}
			set_value(prd, X, 0, tmpmin);
			set_base(prd, X, 0, Ystar);

			tmpmin = INF;
			for (Y = X; Y <= X+K[1]; Y++) {
				tmp = L2(Y, prd);
				if (tmp < tmpmin) {
					tmpmin  = tmp;
					Ystar = Y;
				}
			}
			set_value(prd, X, 1, tmpmin+PES2(X, prd));
			set_base(prd, X, 1, Ystar);

			tmpmin = INF;
			for (Y = X; Y <= X+K[2]; Y++) {
				tmp = L3(Y, prd);
				if (tmp < tmpmin) {
					tmpmin = tmp;
					Ystar = Y;
				}
			}
			set_value(prd, X, 2, tmpmin + PES3(X, prd));
			set_base(prd, X, 2, Ystar);
		}
	}
}

void init()
{
	//TODO: This function initialize all variables including the value
	//	function array.

	int i;
	FILE * fp = fopen("nEchelon.dat", "r");
	fscanf(fp, "%lf%d%d", &beta, &N, &period);
	for (i = 0; i < N; i++) {
		fscanf(fp, "%lf", &h[i]);
	}
	for (i = 0; i < N; i++) {
		fscanf(fp, "%d", &K[i]);
	}
	fscanf(fp, "%lf%d%d%d", &p, &UB, &LB, &D_len);
	for (i = 0; i < D_len; i++) {
		fscanf(fp, "%d", &D[i]);
	}
	for (i = 0; i < D_len; i++) {
		fscanf(fp, "%lf", &P[i]);
	}
}

void Eval_Base()
{
	int X[3], prd, B[3], Y[3];
	for (prd = 1; prd <= period; prd++) {
		for (X[0] = LB; X[0] <= UB; X[0]++) {
			for (X[1] = X[0]; X[1] <= UB; X[1]++) {
				for (X[2] = X[1]; X[2] <= UB; X[2]++) {
					B[0] = get_base(prd, X[0], 0);
					B[1] = get_base(prd, X[1], 1);
					B[2] = get_base(prd, X[2], 2);
					Y[0] = MIN(B[0], X[0]+K[0]);
					Y[0] = MIN(Y[0], X[1]);
					Y[0] = MAX(X[0], Y[0]);
					Y[1] = MIN(B[1], X[1]+K[1]);
					Y[1] = MIN(Y[1], X[2]);
					Y[1] = MAX(X[1], Y[1]);
					Y[2] = MIN(B[2], X[2]+K[2]);
					Y[2] = MAX(X[2], Y[2]);
					set_EvV(prd, X, Ev_L(Y, prd)); 
					if (prd == period) {
						set_policy(X, Y);
					}
				}
			}
		}
	}
}

int main(int argc, const char *argv[])
{
	int i, j, k, x[3], X[3];

	init();

	DP();

	Eval_Base();

	//TODO: read installation inventory and print the optimal base
	while (scanf("%d%d%d", &x[0], &x[1], &x[2]) != EOF) {
		X[0] = x[0]; 	X[1] = x[0] + x[1];	X[2] = X[1]+x[2];
		Ending[0] = get_policy(X, 0);
		Ending[1] = get_policy(X, 1);
		Ending[2] = get_policy(X, 2);
		printf("%d\t%d\t%d\t", X[0], X[1], X[2]);
		printf("%d\t%d\t%d\t", Ending[0], Ending[1], Ending[2]);
		printf("%lf\n", get_EvV(period, X));
	}
	return 0;
}
