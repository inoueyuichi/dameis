/*
 * =====================================================================================
 *
 *       Filename:  K1INF.c
 *
 *    Description:  This is a two-echelon system where the lower installation has 
 *		    an infinite capacity.
 *
 *        Version:  1.0
 *        Created:  2014年11月25日 19时52分00秒
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

//Value Functions
double V[2][MAX_PERIOD][MAX_X * 2];

//Policy
int Plc[2][MAX_X * 2];

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

double f_n_1(int Y_n_1, int prd)
{
	int i;
	double res = 0,
	       tmpexp = 0;
	for (i = 0; i < D_len; i++) {
		if (Y_n_1 > D[i]) {
			res += (h[0] + h[1]) * (Y_n_1 - D[i]) * P[i];
		}
		else {
			res += p * (D[i] - Y_n_1) * P[i];
		}
		tmpexp += V[0][prd-1][MAX_X+Y_n_1 - D[i]] * P[i];
	}
	res -= h[1] * Y_n_1;
	res += beta * tmpexp;
	return res;
}

double f_n_2(int Y_n_2, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res += V[1][prd-1][MAX_X+Y_n_2 - D[i]] * P[i];
	}
	res = res * beta + h[1] * Y_n_2;
	return res;
}

void DP_1(int X1, int prd)
{
	int Y;
	double tmpres;
	for (Y = X1; Y <= UB; Y++) {
		tmpres = f_n_1(Y, prd);
		if (tmpres < V[0][prd][MAX_X+X1]) {
			V[0][prd][MAX_X+X1] = tmpres;
			if (prd == period) {
				Plc[0][MAX_X+X1] = Y;
			}
		}
		else {
			break;
		}
	}
}

void DP_2(int X2, int prd)
{
	double min1 = INF,
	       min2 = INF,
	       tmpres;
	int Y,
	    YUB = X2 + K[1];
	for (Y = X2; Y <= UB; Y ++) {
		tmpres = f_n_2(Y, prd);
		if (tmpres < min1) {
			min1 = tmpres;
		}
		else {
			break;
		}
	}
	for (Y = X2; Y <= YUB; Y++) {
		tmpres = f_n_2(Y, prd);
		if (tmpres < min2) {
			min2 = tmpres;
		}
		else {
			break;
		}
	}
	V[1][prd][MAX_X+X2] = f_n_1(X2, prd) - min1 + min2;
	if (prd == period) {
		Plc[1][MAX_X+X2] = Y;
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
	for (i = 0; i < MAX_PERIOD; i++) {
		for (j = 0; j < MAX_X *2; j++) {
			V[0][i][j] = V[1][i][j] = ((i==0)?0:INF);
		}
	}
}

int main(int argc, const char *argv[])
{
	int i, j, k,
	    X[2];
	init();

	//TODO: k iterates through periods; i and j through all possible states
	for (k = 1; k <= period; k++) {
		for (i = LB; i <= UB; i++) {
			DP_1(i, k);
			DP_2(i, k);
		}
	}

	//TODO: read installation inventory and print the optimal policy
	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2lf\n", X[0], X[1] - X[0],
				X[0], X[1], Plc[0][X[0]+MAX_X] - X[0], 
				Plc[1][X[1]+MAX_X] - X[1], Plc[0][X[0] + MAX_X], 
				Plc[1][X[1] + MAX_X], 
				V[0][period][X[0]+MAX_X] + V[1][period][X[1]+MAX_X]);
	}
	return 0;
}
