/*
 * =====================================================================================
 *
 *       Filename:  2Lagrange.c
 *
 *    Description:  This program uses Lagrangian relaxation to compute the optimal
 *		    policy.
 *
 *        Version:  1.0
 *        Created:  2014年11月24日 13时37分56秒
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
#define MAX_LAMBDA		20
#define N			2
#define INF			99999999
#define MAX(x,y)		((x>y)?x:y)
#define MIN(x,y)		((x<y)?x:y)
//TODO:	Retrieve() is a macro which gives access to V and Plc
#define Retrieve(V,i,X,Y)	\
	V[i][X+MAX_X][Y+MAX_X]

//Value Function
double V[MAX_PERIOD][MAX_X*2][MAX_X*2];

//Policy
//int Plc[MAX_PERIOD][MAX_X*2][MAX_X*2][2];

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

double Jn(int X[], int Y[], double lambda[], int prd)
{
	int i;
	double res, tmpE = 0;
	res = l(Y[0]);
	res += h[1] * (Y[1] - Y[0]);
	res += lambda[0] * (Y[0] - X[0] - K[0]);
	res += lambda[1] * (Y[1] - X[1] - K[1]);
	for (i = 0; i < D_len; i++) {
		tmpE += P[i] * Retrieve(V, prd - 1, Y[0]-D[i], Y[1]-D[i]);
	}
	return (res + beta * tmpE);
}

double Ln(int X[], double lambda[], int prd)
{
	int Y[2];
	double res = INF, tmp;
	for (Y[0] = X[0]; Y[0] <= X[1]; Y[0]++) {
		for (Y[1] = X[1]; Y[1] <= UB; Y[1]++) {
			tmp = Jn(X, Y, lambda, prd);
			if (tmp < res) {
				res = tmp;
			}
		}
	}
	return res;
}

void DP(int X[], int prd)
{
	double lambda[2], tmp;
	for (lambda[0] = 0; lambda[0] <= MAX_LAMBDA; lambda[0] += 1) {
		for (lambda[1] = 0; lambda[1] <= MAX_LAMBDA; lambda[1] += 1) {
			tmp = Ln(X, lambda, prd);
			if (tmp > Retrieve(V, prd, X[0], X[1])) {
				Retrieve(V, prd, X[0], X[1]) = tmp;
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
			for (k = 0; k < MAX_X * 2; k++) {
				V[i][j][k] = ((i==0)?0:-INF);
			}
		}
	}
}

int main(int argc, const char *argv[])
{
	int i, j, k, X[2];
	init();

	//TODO: k iterates through periods; i and j through all possible states
	for (k = 1; k <= period; k++) {
		for (X[0] = LB; X[0] <= UB; X[0]++) {
			for (X[1] = LB; X[1] <= UB; X[1]++) {
				DP(X, k);
			}
		}
	}

	//TODO: read installation inventory and print the optimal policy
	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		printf("%d\t%d\t%.2lf\n", X[0], X[1], Retrieve(V, period, X[0], X[1]));
	}
	return 0;
}
