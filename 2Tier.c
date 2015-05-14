/*
 * =====================================================================================
 *
 *       Filename:  EvalPol.c
 *
 *    Description:  This program evaluate a policy.
 *
 *        Version:  1.0
 *        Created:  2015年4月22日 8时28分37秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LU Tianshu (), tssslu@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>

#define MAX_PERIOD		100
#define MAX_X			200
#define MAX_D_LENGTH		30
#define N			2
#define INF			99999999
#define MAX(x,y)		(((x)>(y))?(x):(y))
#define MIN(x,y)		(((x)<(y))?(x):(y))

//Value Function
double V[MAX_PERIOD][MAX_X*2][MAX_X*2];

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

//Base stock
int S2;
int S1[2];
int XS[2];

//Lower & Upper Bound of X, Number of Periods
int LB, UB, period;

double get_value(int X[], int prd)
{
	return V[prd][X[0]+MAX_X][X[1]+MAX_X];
}

int get_S1(int X2)
{
	if (X2 <= XS[0]) {
		return S1[0];
	}
	else if (X2 >XS[1]) {
		return S1[1];
	}
	else {
		return S1[0]+((double)(S1[1]-S1[0]))
			/((double)(XS[1]-XS[0]))*(X2-XS[0]);
	}
}

void set_value(int X[], int prd, double val)
{
	V[prd][X[0]+MAX_X][X[1]+MAX_X] = val;
}

double l(int Y1)
{
	int i, tmp;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		tmp = Y1 - D[i];
		if (tmp > 0) {
			res += P[i] * (h[0]+h[1]) * tmp;
		}
		else {
			res -= P[i] * p * tmp;
		}
	}
	return res;
}

void policy(int X[], int Y[])
{
	int S1 = get_S1(X[1]);
	Y[0] = MIN(MIN(S1, X[0]+K[0]), X[1]);
	Y[0] = MAX(X[0], Y[0]);
	Y[1] = MIN(S2, X[1]+K[1]);
	Y[1] = MAX(X[1], Y[1]);
}

double Jn(int Y[], int prd)
{
	int i, tmpX[2];
	double res, tmpEV = 0;
	res = l(Y[0]) + h[1] * (Y[1] - Y[0]);
	for (i = 0; i < D_len; i++) {
		tmpX[0] = Y[0] - D[i];
		tmpX[1] = Y[1] - D[i];
		tmpEV += P[i] * get_value(tmpX, prd-1);
	}
	return res + beta * tmpEV;
}

void Vn(int X[], int prd)
{
	int Y[2];
	double tmp;
	policy(X, Y);
	tmp = Jn(Y, prd);
	set_value(X, prd, tmp);
}

void DP()
{
	int prd, X[2];
	for (prd = 1; prd <= period; prd++) {
		for (X[0] = LB; X[0] <= UB; X[0]++) {
			for (X[1] = LB; X[1] <= UB; X[1]++) {
				Vn(X, prd);
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
}

int main(int argc, const char *argv[])
{
	int X[2];
	init();
	scanf("%d%d%d%d%d", &S1[0], &XS[0], &S1[1], &XS[1], &S2);
	DP();
	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		printf("%d\t%d\t%.2lf\n", X[0], X[1],
				get_value(X, period));
	}
	return 0;
}
