/*
 * =====================================================================================
 *
 *       Filename:  9Targets.c
 *
 *    Description:  This program gets the nine target levels described in Ji 2015
 *
 *        Version:  1.0
 *        Created:  2015年05月14日 12时49分24秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LU Tianshu (), tssslu@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <assert.h>
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
double V[MAX_PERIOD][MAX_X*2][MAX_X*2];

//Policy
int Plc[MAX_X*2][MAX_X*2][2];

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

double get_value(int X1, int X2, int t)
{
	return V[t][X1+MAX_X][X2+MAX_X];
}

void set_value(int X1, int X2, int t, double val)
{
	V[t][X1+MAX_X][X2+MAX_X] = val;
}

int get_policy(int X1, int X2, int Y[])
{
	Y[0] = Plc[X1+MAX_X][X2+MAX_X][0];
	Y[1] = Plc[X1+MAX_X][X2+MAX_X][1];
}

void set_policy(int X1, int X2, int Y[])
{
	Plc[X1+MAX_X][X2+MAX_X][0] = Y[0];
	Plc[X1+MAX_X][X2+MAX_X][1] = Y[1];
}

double phi(int Y1)
{
	int i;
	double tmp = 0, res = 0;
	for (i = 0; i < D_len; i++) {
		if (Y1 - D[i] >= 0) {
			tmp = (h[0]+h[1]) * (Y1-D[i]);
		}
		else {
			tmp = p * (D[i] - Y1);
		}
		res += P[i] * tmp;
	}
	return res;
}

double L(int X[], int Y[])
{
	return phi(Y[0]) + h[1] * (Y[1]-Y[0]);
}

double J(int X[], int Y[], int t)
{
	int i;
	double exp = 0;
	for (i = 0; i < D_len; i++) {
		exp+=P[i] * get_value(X[0]-D[i], X[1]-D[i], t);
	}
	return L(X, Y) + beta * exp;
}

void DP()
{
	int t, X[2], Y[2], tmpJ;
	for (t = 1; t < period; t++) {
		for (X[0] = LB; X[0] < UB; X[0]++) {
			for (X[1] = LB; X[1] < UB; X[1]++) {
				tmpJ = INF;
				for (Y[0] = 0; Y[0] < MIN(X[1], X[0]+K[0]); Y[0]++) {
					for (Y[1] = 0; Y[1] < X[1]+K[1]; Y[1]++) {
						if (J(X,Y,t) < tmpJ) {
							tmpJ = J(X,Y,t);
						}
					}
				}
				set_value(X[0], X[1], t, tmpJ);
				if (t == period) {
					set_policy(X[0], X[1], Y);
				}
			}
		}
	}
}

int Y11(int X[], int t)
{
	int Y[2], res;
	double tmp = INF;
	Y[1] = X[1];
	for (Y[0] = LB; Y[0] < UB; Y[0]++) {
		if (J(Y,X,t) < tmp) {
			tmp = J(Y, X, t);
			res = Y[0];
		}
	}
	return res;
}

int Y12(int X[], int t)
{
	int Y[2], res;
	double tmp = INF;
	for (Y[0] = LB; Y[0] < UB; Y[0]++) {
		Y[1] = Y[0] + (t-1)*K[0] -K[1];
		if (J(Y,X,t) < tmp) {
			tmp = J(Y,X,t);
			res = Y[0];
		}
	}
	return res;
}

int Y13(int X[], int t)
{
	int Y[2], res;
	double tmp = INF;
	for (Y[0] = LB; Y[0] < UB; Y[0]++) {
		Y[1] = Y[0] + (t-2) * K[0];
		if (J(Y,X,t) < tmp) {
			tmp = J(Y,X,t);
			res = Y[0];
		}
	}
	return res;
}

int Y14(int X[], int t)
{
	int Y[2], res;
	double tmp = INF;
	for (Y[0] = LB; Y[0] < UB; Y[0]++) {
		for (Y[1] = LB; Y[1] < UB; Y[1]++) {
			if (J(Y,X,t) < tmp) {
				tmp = J(Y, X, t);
				res = Y[0];
			}
		}
	}
	return res;
}

int Y15(int X[], int t)
{
	int Y[2], res;
	double tmp = INF;
	for (Y[0] = LB; Y[0] < UB; Y[0]++) {
		Y[1] = Y[0] + (t-2)*K[0]-K[1];
		if (J(Y,X,t) < tmp) {
			tmp = J(Y,X,t);
			res = Y[0];
		}
	}
	return res;
}

int Y21(int X[], int t)
{
	return Y12(X, t) +(t-1)*K[0]-K[1];
}

int Y22(int X[], int t)
{
	return Y13(X, t) + (t-2)*K[0];
}

int Y23(int X[], int t)
{
	int Y[2], res;
	double tmp = INF;
	for (Y[0] = LB; Y[0] < UB; Y[0]++) {
		for (Y[1] = LB; Y[1] < UB; Y[1]++) {
			if (J(Y,X,t) < tmp) {
				tmp = J(Y, X, t);
				res = Y[1];
			}
		}
	}
	return res;
}

int Y24(int X[], int t)
{
	return Y15(X, t)+(t-2)*K[0]-K[1];
}

void find_intersect(int X[], int t[])
{
	int i, Y[2];
	for (i = 0; i < period; i++) {
		Y[0] = X[1]+K[1]-(i-1)*K[0];
		if (X[0] <= Y[0] && Y[0] <= X[0]+K[0]) {
			break;
		}
	}
	if (i ==  period+1) {
		t[0] = t[1] = t[2] = INF;
		return;
	}
	t[0] = i;
	Y[0] = X[1]-(i-2)*K[0];
	if (X[0] <= Y[0] && Y[0] <=X[0]+K[0]) {
		t[1] = i;
	}
	else {
		t[1] = t[2] = INF;
		return;
	}
	Y[0] = X[1] +K[1]-(i-2)*K[0];
	if (X[0] <= Y[0] && Y[0] <=X[0]+K[0]) {
		t[2] = i;
	}
	else {
		t[2] = INF;
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
	for (i = 0; i < MAX_PERIOD; i++) {
		for (j = 0; j < MAX_X *2; j++) {
			for (k = 0; k < MAX_X * 2; k++) {
				V[i][j][k] = ((i==0)?0:INF);
			}
		}
	}
}

int main(int argc, const char *argv[])
{
	int X[2], Y[2], t[2];
	init();

	DP();

	while (scanf("%d%d", &X[0], &X[1])!=EOF) {
		X[1] += X[0];
		get_policy(X[0], X[1], Y);
		printf("%d\t%d\t%d\t%d\t", X[0], X[1], Y[0], Y[1]);
		printf("%.2lf\t", get_value(X[0], X[1], period));
		printf("%d\t", Y11(X, period));
		find_intersect(X, t);
		if (t[0] != INF) {
			printf("%d\t", Y12(X, t[0]));
		}
		else {
			printf("NA\t");
		}
		if (t[1] != INF) {
			printf("%d\t", Y13(X, t[1]));
		}
		else {
			printf("NA\t");
		}
		printf("%d\t", Y14(X, period));
		if (t[2] != INF) {
			printf("%d\t", Y15(X, t[2]));
		}
		else {
			printf("NA\t");
		}
		if (t[0] != INF) {
			printf("%d\t", Y21(X, t[0]));
		}
		else {
			printf("NA\t");
		}
		if (t[1] != INF) {
			printf("%d\t", Y22(X, t[1]));
		}
		else {
			printf("NA\t");
		}
		printf("%d\t", Y23(X, period));
		if (t[2] != INF) {
			printf("%d\n", Y24(X, t[2]));
		}
		else {
			printf("NA\n");
		}
	}
	return 0;
}
