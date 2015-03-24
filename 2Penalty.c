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

#define MAX_PERIOD		100
#define MAX_X			200
#define MAX_D_LENGTH		30
#define N			2
#define INF			99999999
#define MAX(x,y)		((x>y)?x:y)
#define MIN(x,y)		((x<y)?x:y)

//Value Function
double V[MAX_PERIOD][MAX_X*2][2];

//Policy
int Plc[MAX_X*2][2];

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

//Capacity Constraint
int K[N];

//Discount Factor
double beta;

//Lower & Upper Bound of X, Number of Periods
int LB, UB, period;

//Initial Installation & Echelon Inventory
int x[2], X[2];

void set_value(int prd, int X, int cnt, double val)
{
	V[prd][X+MAX_X][cnt] = val;
}

double get_value(int prd, int X, int cnt)
{
	return V[prd][X+MAX_X][cnt];
}

void set_policy(int X, int cnt, int val)
{
	Plc[X+MAX_X][cnt] = val;
}

int get_policy(int X, int cnt)
{
	return Plc[X+MAX_X][cnt];
}

double l(int Y1)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		if (Y1 > D[i]) {
			res += P[i] * (h[0]+h[1]) * (Y1-D[i]);
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
	return l(Y1) - h[1]*Y1 + beta*EV;
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

double PCC1(int X1, int prd)
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
	if (Ystar <= X1 + K[0]) {
		return 0;
	}
	else {
		return L1(X1+K[0], prd) - tmpmin;
	}
}

double PCC2(int X2, int prd)
{
	int Ystar, Y2;
	double tmp, tmpmin = INF;
	for (Y2 = LB; Y2 <= UB; Y2++) {
		tmp = L2(Y2, prd);
		if (tmp < tmpmin) {
			tmpmin = tmp;
			Ystar = Y2;
		}
	}
	if (Ystar <= X2 + K[1]) {
		return 0;
	}
	else {
		return L2(X2+K[1], prd) - tmpmin;
	}
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

void DP()
{
	int prd, X, Y, Ystar;
	double tmpV, tmpmin, tmp;
	for (prd = 1; prd <= period; prd++) {
		for (X = LB; X <= UB; X++) {
			tmpV = PCC1(X, prd);
			tmpmin = INF;
			for (Y = X; Y <= UB; Y++) {
				tmp = L1(Y, prd);
				if (tmp < tmpmin) {
					tmpmin = tmp;
					Ystar = Y;
				}
			}
			set_value(prd, X, 0, tmpmin+tmpV);
			if (prd == period) {
				set_policy(X, 0, Ystar);
			}

			tmpV = PCC2(X, prd) + PES2(X, prd);
			tmpmin = INF;
			for (Y = X; Y <= UB; Y++) {
				tmp = L2(Y, prd);
				if (tmp < tmpmin) {
					tmpmin  = tmp;
					Ystar = Y;
				}
			}
			set_value(prd, X, 1, tmpmin+tmpV);
			if (prd == period) {
				set_policy(X, 1, Ystar);
			}
		}
	}
}

void init()
{
	//TODO: This function initialize all variables including the value
	//	function array.
	//	This part can be changed in order to read data from a file.

	int i;
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
	for (i = 0; i < MAX_PERIOD*2; i++) {
		V[0][i][0] = V[0][i][1] = 0;
	}
}

int main(int argc, const char *argv[])
{
	int i, j, k;
	init();

	DP();

	//TODO: read installation inventory and print the optimal policy
	while (scanf("%d%d", &x[0], &x[1]) != EOF) {
		X[0] = x[0]; 	X[1] = x[0] + x[1];
		Ending[0] = get_policy(X[0], 0);
		Ending[1] = get_policy(X[1], 1);
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", x[0], x[1],
			X[0], X[1], Ending[0] - X[0], Ending[1] - X[1],
			Ending[0], Ending[1], 
			get_value(period, X[0], 0)+get_value(period,X[1], 1));
	}
	return 0;
}
