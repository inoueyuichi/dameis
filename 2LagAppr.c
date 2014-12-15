/*
 * =====================================================================================
 *
 *       Filename:  2LagAppr.c
 *
 *    Description:  This program implements the approximate algorithm of a 2-echelon
 *			inventory system.
 *
 *        Version:  1.0
 *        Created:  2014年12月09日 08时53分48秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LU Tianshu (), tssslu@gmail.com
 *   Organization:  THU IE
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

//Value Functions
double V[N][MAX_PERIOD][MAX_X*2];

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

double get_value(int ech, int prd, int X)
{
	return V[ech][prd][X+MAX_X];
}

void set_value(int ech, int prd, int X, double uv)
{
	V[ech][prd][X+MAX_X] = uv;
}

double l(int Y1)
{
	int i;
	double res = 0, tmp;
	for (i = 0; i < D_len; i++) {
		if (Y1 > D[i]) {
			tmp = (h[0]+h[1]) * (Y1-D[i]);
		}
		else {
			tmp = p * (D[i] - Y1);
		}
		res += P[i] * tmp;
	}
	return res;
}

double J1(int Y1, int prd)
{
	int i;
	double tmpEV = 0;
	for (i = 0; i < D_len; i++) {
		tmpEV += get_value(0, prd-1, Y1-D[i]) * P[i];
	}
	return l(Y1) - h[1] * Y1 + beta * tmpEV;
}

void DP1(int X1, int prd)
{
	int i, Y1, Y1U;
	double tmpJ, tmpmin = INF;
	Y1U = X1 + K[0];
	for (Y1 = X1; Y1 <= Y1U; Y1 ++) {
		tmpJ = J1(Y1, prd);
		if (tmpJ < tmpmin) {
			tmpmin = tmpJ;
		}
		else {
			break;
		}
	}
	set_value(0, prd, X1, tmpmin);
}

double J2(int Y2, int prd)
{
	int i;
	double tmpEV = 0;
	for (i = 0; i < D_len; i++) {
		tmpEV += get_value(1, prd-1, Y2 - D[i]) * P[i];
	}
	return h[1] * Y2 + beta * tmpEV;
}

void DP2(int X2, int prd) 
{
	int i, Y1, Y2, Y2U;
	double res, tmpEV = 0, tmpJ,
	       tmpmin = INF;
	Y2U = X2 + K[1];
	res = J1(X2, prd);
	for (Y2 = X2; Y2 <= Y2U; Y2 ++) {
		tmpJ = J2(Y2, prd);
		if (tmpJ < tmpmin) {
			tmpmin = tmpJ;
		}
		else {
			break;
		}
	}
	res += tmpmin;
	tmpmin = INF;
	for (Y1 = X2; 1; Y1 ++) {
		tmpJ = J1(Y1, prd);
		if (tmpJ < tmpmin) {
			tmpmin = tmpJ;
		}
		else {
			break;
		}
	}
	set_value(1, prd, X2, res - tmpmin);
}

void get_policy(int X[], int Y[])
{
	int Y1U, Y2U, Y1, Y2;
	double tmpJ, tmpmin = INF;
	Y1U = MIN(X[1], X[0] + K[0]);
	Y2U = X[1] + K[1];
	for (Y1 = X[0]; Y1 <= Y1U; Y1++) {
		tmpJ = J1(Y1, period);
		if (tmpJ < tmpmin) {
			tmpmin = tmpJ;
			Y[0] = Y1;
		}
		else {
			break;
		}
	}
	tmpmin = INF;
	for (Y2 = X[1]; Y2 <= Y2U; Y2 ++) {
		tmpJ = J2(Y2, period);
		if (tmpJ < tmpmin) {
			tmpmin = tmpJ;
			Y[1] = Y2;
		}
		else {
			break;
		}
	}
}

void init()
{
	//TODO: This function initialize all variables including the value
	//	function array.
	//	This part can be changed in order to read data from a file.

	int i, j;
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
	for (i = 0; i < MAX_X * 2; i++) {
		V[0][0][i] = V[1][0][i] = 0;
	}
	for (i = 0; i <= period; i++) {
		for (j = 0; j < 2*MAX_X; j++) {
			V[0][i][j] = V[1][i][j] = (i==0)?0:INF;
		}
	}
}

int main(int argc, const char *argv[])
{
	int i, j, X[2], Y[2];
	init();
	for (i = 1; i <= period; i++) {
		for (j = LB; j <= UB; j++) {
			DP1(j, i);
			DP2(j, i);
		}
	}

	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		get_policy(X, Y);
		printf("%d\t%d\t%d\t%d\n", X[0], X[1], Y[0], Y[1]);
	}
	return 0;
}
