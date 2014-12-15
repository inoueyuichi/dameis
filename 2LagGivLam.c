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
double V[MAX_PERIOD][MAX_X*2][MAX_X*2];

//Policy
int Plc[MAX_X*2][MAX_X*2][2];

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
double lambda[2];

double get_value(int X[], int prd)
{
	return V[prd][X[0]+MAX_X][X[1]+MAX_X];
}

void set_value(int X[], int prd, int val)
{
	V[prd][X[0]+MAX_X][X[1]+MAX_X] = val;
}

void set_policy(int X[], int Y[])
{
	Plc[X[0]+MAX_X][X[1]+MAX_X][0] = Y[0];
	Plc[X[0]+MAX_X][X[1]+MAX_X][1] = Y[1];
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

double Jn(int X[], int Y[], int prd)
{
	int i, tmpX[2];
	double res, tmpEV = 0;
	res = l(Y[0]) + h[1] * (Y[1] - Y[0]) + 
		lambda[0] * (Y[0] - X[0] - K[0]) +
		lambda[1] * (Y[1]- X[1] - K[1]);
	for (i = 0; i < D_len; i++) {
		tmpX[0] = Y[0] - D[i];
		tmpX[1] = Y[1] - D[i];
		tmpEV += P[i] * get_value(tmpX, prd - 1);
	}
	return res + beta * tmpEV;
}

double Ln(int X[], int prd)
{
	int Y[2];
	double tmpJ, tmpmin = INF;
	for (Y[0] = X[0]; Y[0] <= X[1]; Y[0] ++) {
		for (Y[1] = X[1]; Y[1] <= UB; Y[1] ++) {
			tmpJ = Jn(X, Y, prd);
			if (tmpJ < tmpmin) {
				tmpmin = tmpJ;
				if (prd == period) {
					set_policy(X, Y);
				}
			}
		}
	}
	return tmpmin;
}

void DP()
{
	int prd, X[2];
	for (prd = 1; prd <= period; prd++) {
		for (X[0] = LB; X[0] <= UB; X[0] ++) {
			for (X[1] = LB; X[1] <= UB; X[1] ++) {
				set_value(X, prd, 
					Ln(X, prd));
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
				V[i][j][k] = ((i==0)?0:INF);
			}
		}
	}
}

int main(int argc, const char *argv[])
{
	int X[2];
	init();
	scanf("%lf%lf", &lambda[0], &lambda[1]);
	DP();
	while (scanf("%d%d", &X[0], &X[1]) != EOF) {
		X[1] += X[0];
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2lf\n", X[0], X[1]-X[0], X[0], X[1], 
				Plc[X[0]+MAX_X][X[1]+MAX_X][0],
				Plc[X[0]+MAX_X][X[1]+MAX_X][1], 
				Plc[X[0]+MAX_X][X[1]+MAX_X][0]-X[0],
				Plc[X[0]+MAX_X][X[1]+MAX_X][1]-X[1],
				get_value(X, period));
	}
	return 0;
}
