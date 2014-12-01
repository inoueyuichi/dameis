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
#define MIN(x,y)		\
	((x<y)?x:y)
#define MAX(x,y)		\
	((x>y)?x:y)

//Value Function
double V[2][MAX_PERIOD][MAX_X*2];
int Plc[2][MAX_X*2];
#define V(ech,prd,x)		\
	(V[ech-1][prd][MAX_X+x])
#define Plc(ech,x)		\
	(Plc[ech-1][MAX_X+x])

//Incremental Holding Cost
double h[N];

//Penalty Cost
double p;

//Demand Distribution and Expected Demand
int D_len;
int D[MAX_D_LENGTH];
double P[MAX_D_LENGTH];

//Discount Factor
double beta;

//Lower & Upper Bound of X, Number of Periods
int LB, UB, period;

double l(int x)
{
	int i, H = h[0]+h[1];
	double res = 0, tmp;
	for (i = 0; i < D_len; i++) {
		tmp = x - D[i];
		if (tmp > 0) {
			res += P[i] * H * tmp;
		}
		else {
			res -= P[i] * p * tmp;
		}
	}
	return res;
}

double delta(int X2, int X1_star, int prd)
{
	int i;
	double res, tmp = 0;
	if (X2 >= X1_star) {
		return 0;
	}
	res = l(X2) -h[1]*(X2-X1_star) - l(X1_star);
	for (i = 0; i < D_len; i++) {
		tmp += P[i] * V(1, prd-1, X2-D[i]);
		tmp -= P[i] * V(1, prd-1, X1_star-D[i]);
	}
	return res + beta * tmp;
}

double J1(int Y1, int prd)
{
	int i;
	double tmp = 0;
	for (i = 0; i < D_len; i++) {
		tmp += P[i] * V(1, prd-1, Y1-D[i]);
	}
	return l(Y1) - h[1] * Y1 + beta * tmp;
}

double J2(int Y2, int prd)
{
	int i;
	double tmp = 0;
	for (i = 0; i < D_len; i++) {
		tmp += P[i] * V(2,prd-1, Y2-D[i]);
	}
	return h[1]*Y2 + beta * tmp;
}

int DP1(int X1, int prd)
{
	int Y1, tmpPlc;
	double tmpJ;
	for (Y1 = X1; 1; Y1 ++) {
		tmpJ = J1(Y1, prd);
		if (tmpJ < V(1, prd, X1)) {
			V(1, prd, X1) = tmpJ;
			tmpPlc = Y1;
		}
		else {
			break;
		}
	}
	return tmpPlc;
}

int DP2(int X2, int X1_star, int prd)
{
	int Y2, tmpPlc;
	double tmpJ;
	for (Y2 = X2; 1; Y2 ++) {
		tmpJ = J2(Y2, prd);
		tmpJ += delta(X2, X1_star, prd);
		if (tmpJ < V(2, prd, X2)) {
			V(2, prd, X2) = tmpJ;
			tmpPlc = Y2;
		}
		else {
			break;
		}
	}
	return tmpPlc;
}

void DP(int prd)
{
	int X1, X2, tmpPlc, X1_star = INF;
	for (X1 = LB; X1 <= UB; X1++) {
		tmpPlc = DP1(X1, prd);
		if (tmpPlc < X1_star) {
			X1_star = tmpPlc;
		}
		if (prd == period) {
			Plc(1, X1) = tmpPlc;
		}
	}
	for (X2 = LB; X2 <= UB; X2++) {
		tmpPlc = DP2(X2, X1_star, prd);
		if (prd == period) {
			Plc(2, X2) = tmpPlc;
		}
	}
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
	for (i = 0; i < N; i++) {
		for (j = 0; j < MAX_PERIOD; j++) {
			for (k = 0; k < MAX_X*2; k++) {
				V[i][j][k] = (j == 0) ? 0: INF;
			}
		}
	}
}

int main(int argc, const char *argv[])
{
	int prd, X[2], E[2];
	init();

	for (prd = 1; prd <= period; prd++) {
		DP(prd);
	}

	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		E[0] = MAX(X[0], MIN(X[1], Plc(1, X[0])));
		E[1] = MAX(X[1], Plc(2, X[1]));
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", X[0], X[1] - X[0],
				X[0], X[1], E[0] - X[0], E[1] - X[1], E[0], E[1],
				V(1, period, X[0])+V(2,period, X[1]));
	}
	return 0;
}
