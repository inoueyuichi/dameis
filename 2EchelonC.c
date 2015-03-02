/*
 * =====================================================================================
 *
 *       Filename:  2Echelon.c
 *
 *    Description:  Model in Parker and Kapuscinski 2004, page 741.
 *
 *        Version:  1.0
 *        Created:  2014年11月07日 13时38分17秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  LU Tianshu (), tssslu@gmail.com
 *   Organization:  THU IE
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
//TODO:	Retrieve() is a macro which gives access to V and Plc
#define Retrieve(V,i,X,Y)	\
	V[i][X+MAX_X][Y+MAX_X]
#define ASSERT(Y)		\
	(assert(Y[0] >= -MAX_X && Y[0] <= MAX_X &&\
		Y[1] >= -MAX_X && Y[1] <= MAX_X))

//Value Function
double V[MAX_PERIOD][MAX_X*2][MAX_X*2];

//Policy
int Plc[MAX_PERIOD][MAX_X*2][MAX_X*2][2];

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

//Initial Installation & Echelon Inventory
int x[2], X[2];

double L(int Y[])
{
	//TODO: Given a decison Y[N], this function tells its holding 
	//	and penalty costs.
	//VAR:	EDYp		E[max(D-Y,0)]
	//	pen_hol		p+h_1+h_2+...+h_N

	int i;
	double res = 0, pen_hol = p, EDYp = 0;
	for (i = 0; i < D_len; i++) {
		EDYp += MAX(D[i]-Y[0], 0) * P[i];
	}
	for (i = 0; i < N; i++) {
		res += h[i] * (Y[i] - ED);
		pen_hol += h[i];
	}
	res += pen_hol * EDYp;
	return res;
}

double J(int Y[], int prd)
{
	//TODO: Given a decision Y[N] and the current period, this
	//	function tells the value of objective function.
	//VAR:	EV		E[V_n-1(Y-D)]

	int i;
	double res, EV = 0;
	ASSERT(Y);
	res = L(Y);
	for (i = 0; i < D_len; i++) {
		EV += P[i] * 
			(Retrieve(V, prd - 1, Y[0] - D[i], Y[1] - D[i]));
	}
	res += beta * EV;
	return res;
}

void DP(int X0, int X1, int prd)
{
	//TODO:	Given the current period prd, the state variables X0 and X1,
	//	this function iterates the value function for one step.
	//VAR:	Y0U		Upper Bound of Y0 = min(X1, X0+K_0)
	//	Y1U		Upper Bound of Y1 = X1+K_1

	int  Y1L, Y1U, Y[2];
	double tmpJ;
	Y1U = X1 + K[1];
	for (Y[0] = X0; Y[0] <= X1; Y[0] ++) {
		Y1L = MAX(X1, Y[0]+K[0]-K[1]);
		Y1U = MIN(X[1]+K[1], Y0+K[0]);
		for (Y[1] = Y1L; Y[1] <= Y1U; Y[1] ++) {
			tmpJ = J(Y, prd);
			if (Retrieve(V, prd, X0, X1) > tmpJ) {
				Retrieve(V, prd, X0, X1) = tmpJ;
				Retrieve(Plc, prd, X0, X1)[0] = Y[0];
				Retrieve(Plc, prd, X0, X1)[1] = Y[1];
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
	int i, j, k;
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
	while (scanf("%d%d", &x[0], &x[1]) != EOF) {
		X[0] = x[0]; 	X[1] = x[0] + x[1];
		Ending[0] = Retrieve(Plc, period, X[0], X[1])[0];
		Ending[1] = Retrieve(Plc, period, X[0], X[1])[1];
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\n", x[0], x[1],
				X[0], X[1], Ending[0] - X[0], Ending[1] - X[1],
				Ending[0], Ending[1], Retrieve(V, period, X[0], X[1]));
	}
	return 0;
}
