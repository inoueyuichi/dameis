/*
 * =====================================================================================
 *
 *       Filename:  InfCapa.c
 *
 *    Description:  This program implements the algorithm in Clark & Scarf, 1960
 *
 *        Version:  1.0
 *        Created:  2014年11月27日 18时13分43秒
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

//Value Function
double V[2][MAX_PERIOD][MAX_X*2];

//Incremental Holding Cost
double h[N];

//Penalty Cost
double p;

//Demand Distribution and Expected Demand
int D_len;
int D[MAX_D_LENGTH];
double P[MAX_D_LENGTH];
double ED = 0;

//Discount Factor
double beta;

//Lower & Upper Bound of X, Number of Periods
int LB, UB, period;

double Jn1(int Yn1, int prd)
{
	int i;
	double res, tmpexp = 0, tmp;
	res = h[0] * (Yn1 + ED);
	for (i = 0; i < D_len; i++) {
		tmp = D[i] - Yn1;
		tmpexp += (tmp > 0) ? P[i]*tmp : 0;
	}
	res += (p+h[0]+h[1]) * tmpexp;
	tmpexp = 0;
	for (i = 0; i < D_len; i++) {
		tmpexp += P[i]*V[0][prd-1][Yn1-D[i]+MAX_X];
	}
	return res + beta * tmpexp;
}

double Jn2(int Yn2, int prd)
{
	int i;
	double res, tmpexp = 0;
	res = h[1] * (Yn2 - ED);
	for (i = 0; i < D_len; i++) {
		tmpexp += P[i]*V[1][prd-1][Yn2-D[i]+MAX_X];
	}
	return res + beta * tmpexp;
}

int policy1(int Xn1, int Xn2, int prd)
{
	int zn1 = Xn2, Yn1;
	double tmp, tmpJ = INF;
	for (Yn1 = Xn1; Yn1 <= Xn2; Yn1++) {
		tmp = Jn1(Yn1, prd);
		if (tmp < tmpJ) {
			tmpJ = tmp;
			zn1 = Yn1;
		}
		else {
			break;
		}
	}
	return (zn1 < Xn2) ? zn1 : Xn2;
}

int policy2(int Xn2, int prd)
{
	int zn2, Yn2;
	double tmp, tmpJ = INF;
	for (Yn2 = Xn2; 1; Yn2++) {
		tmp = Jn2(Yn2, prd);
		if (tmp < tmpJ) {
			tmpJ = tmp;
			zn2 = Yn2;
		}
		else {
			break;
		}
	}
	return zn2;
}

void DP(int Xn1, int Xn2, int prd)
{
	int Yn1, Yn2;
	Yn1 = policy1(Xn1, Xn2, prd);
	Yn2 = policy2(Xn2, prd);
	V[0][prd][Xn1+MAX_X] = Jn1(Yn1, prd);
	V[1][prd][Xn2+MAX_X] = Jn2(Yn2, prd);
}

void init()
{
	//TODO: This function initialize all variables including the value
	//	function array.
	//	This part can be changed in order to read data from a file.

	int i, j, k, K[2];
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
			V[0][i][j] = V[1][i][j] = 
				(i == 0) ? 0 : INF;
		}
	}
}

int main(int argc, const char *argv[])
{
	int i, j, k, X[2], Ending[2];
	double cost;
	init();

	//TODO: k iterates through periods; i and j through all possible states
	for (k = 1; k <= period; k++) {
		for (i = LB; i <= UB; i++) {
			for (j = LB; j <= UB; j++) {
				DP(i, j, k);
			}
		}
	}

	//TODO: read installation inventory and print the optimal policy
	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		Ending[0] = policy1(X[0], X[1], period);
		Ending[1] = policy2(X[1], period);
		cost = V[0][period][X[0]+MAX_X] + V[1][period][X[1]+MAX_X];
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", X[0], X[1]-X[0],
				X[0], X[1], Ending[0] - X[0], Ending[1] - X[1],
				Ending[0], Ending[1], cost);
	}
	return 0;
}
