/*
 * =====================================================================================
 *
 *       Filename:  LowerBounds.c
 *
 *    Description:  This program calculates the lower bound of a three-echelon
 *			inventory system.
 *
 *        Version:  1.0
 *        Created:  2015年06月03日 20时25分45秒
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
#define N			3
#define INF			99999999
#define MAX(x,y)		((x>y)?x:y)
#define MIN(x,y)		((x<y)?x:y)

//Value Function
double U1[3][MAX_PERIOD][MAX_X*2];
double U2[3][MAX_PERIOD][MAX_X*2];

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

double get_U11(int X1, int prd)
{
	return U1[0][prd][X1+MAX_X];
}

double get_U21(int X1, int prd)
{
	return U2[0][prd][X1+MAX_X];
}

void set_U11(int X1, int prd, double val)
{
	U1[0][prd][X1+MAX_X] = val;
}

void set_U21(int X1, int prd, double val)
{
	U2[0][prd][X1+MAX_X] = val;
}

double get_U12(int X2, int prd)
{
	return U1[1][prd][X2+MAX_X];
}

double get_U22(int X2, int prd)
{
	return U2[1][prd][X2+MAX_X];
}

void set_U12(int X2, int prd, double val)
{
	U1[1][prd][X2+MAX_X] = val;
}

void set_U22(int X2, int prd,int val)
{
	U2[1][prd][X2+MAX_X] = val;
}

double get_U13(int X3, int prd)
{
	return U1[2][prd][X3+MAX_X];
}

double get_U23(int X3, int prd)
{
	return U2[2][prd][X3+MAX_X];
}

void set_U13(int X3, int prd, double val)
{
	U1[2][prd][X3+MAX_X] = val;
}

void set_U23(int X3, int prd, double val)
{
	U2[2][prd][X3+MAX_X] = val;
}

double L1(int Y1)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res += (h[0]+h[1]+h[2])*MAX(Y1-D[i],0);
		res += p * MAX(D[i]-Y1, 0);
	}
	return res -(h[1]+h[2])*Y1;
}

double L2(int Y2)
{
	return h[1]*Y2;
}

double L3(int Y3)
{
	return h[2]*Y3;
}

double J11(int Y1, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res+=get_U11(Y1 - D[i], prd);
	}
	return L1(Y1)+res*beta;
}

double J21(int Y1, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res += get_U21(Y1-D[i], prd);
	}
	return L1(Y1)+res*beta;
}

double J12(int Y2, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res+=get_U12(Y2-D[i], prd);
	}
	return L2(Y2) + res*beta;
}

double J22(int Y2, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res += get_U22(Y2-D[i], prd);
	}
	return L2(Y2)+res*beta;
}

double J13(int Y3, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res+=get_U13(Y3-D[i], prd);
	}
	return L3(Y3)+res*beta;
}

double J23(int Y3, int prd)
{
	int i;
	double res = 0;
	for (i = 0; i < D_len; i++) {
		res += get_U23(Y3-D[i], prd);
	}
	return L3(Y3)+beta*res;
}

void update_U11(int prd)
{
	int X1, Y1;
	double tmpmin, tmp;
	for (X1 = LB; X1 <= UB; X1++) {
		tmpmin = INF;
		for (Y1 = X1; Y1 <= UB; Y1++) {
			tmp = J11(Y1, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
			}
		}
		set_U11(X1, prd, tmpmin);
	}
}

void update_U12(int prd)
{
	int X2, Y2, i;
	double tmpmin, tmp;
	for (X2 = LB; X2 <= UB; X2++) {
		tmpmin = INF;
		for (Y2 = X2; Y2 < UB; Y2++) {
			tmp = J12(Y2, prd-1)-J11(Y2, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
				tmp += L1(X2);
				for (i = 0; i < D_len; i++) {
					tmp += beta*get_U11(X2-D[i], prd-1);
				}
			}
		}
		set_U12(X2, prd, tmpmin);
	}
}

void update_U13(int prd)
{
	int X3, Y3, i, Y2;
	double tmpmin, tmp, res;
	for (X3 = LB; X3 <= UB; X3++) {
		tmpmin = INF;
		for (Y3 = X3; Y3 <= MIN(UB, X3+K[2]); Y3++) {
			tmp = J13(Y3, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
			}
		}
		res = tmpmin;
		tmpmin = INF;
		for (Y2 = X3; Y2 <= UB; Y2++) {
			tmp = J12(Y2, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
			}
		}
		res -= tmpmin;
		res += L2(X3);
		for (i = 0; i < D_len; i++) {
			res += beta * get_U12(X3-D[i], prd-1);
		}
		set_U13(X3, prd, res);
	}
}

void DP_U1()
{
	int prd = 0;
	for (prd = 1; prd <= period; prd++) {
		update_U11(prd);
		update_U12(prd);
		update_U13(prd);
	}
}

double get_U1(int X[])
{
	int res;
	res = get_U11(X[0], period);
	res += get_U12(X[1], period);
	res += get_U13(X[2], period);
	return res;
}

void update_U21(int prd)
{
	int X1, Y1;
	double tmpmin, tmp;
	for (X1 = LB; X1 <= UB; X1++) {
		tmpmin = INF;
		for (Y1 = X1; Y1 <= MIN(UB, X1+K[0]); Y1++) {
			tmp = J21(Y1, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
			}
		}
		set_U21(X1, prd, tmpmin);
	}
}

void update_U22(int prd)
{
	int X2, Y2;
	double tmpmin, tmp;
	for (X2 = LB; X2 <= UB; X2++) {
		tmpmin = INF;
		for (Y2 = X2; Y2 < MIN(UB, X2+K[1]); Y2++) {
			tmp = J22(Y2, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
			}
		}
		set_U22(X2, prd, tmpmin);
	}
}

void update_U23(int prd)
{
	int X3, Y3;
	double tmpmin, tmp;
	for (X3 = LB; X3 <= UB; X3++) {
		tmpmin = INF;
		for (Y3 = X3; Y3 < MIN(UB, X3+K[2]); Y3++) {
			tmp = J23(Y3, prd-1);
			if (tmp < tmpmin) {
				tmpmin = tmp;
			}
		}
		set_U23(X3, prd, tmpmin);
	}
}

void DP_U2()
{
	int prd = 0;
	for (prd = 1; prd <= period; prd++) {
		update_U21(prd);
		update_U22(prd);
		update_U23(prd);
	}
}

double get_U2(int X[])
{
	int res;
	res = get_U21(X[0], period);
	res += get_U22(X[1], period);
	res += get_U23(X[2], period);
	return res;
}

void init()
{
	//TODO: This function initialize all variables including the value
	//	function array.

	int i;
	FILE * fp = fopen("nEchelon.dat", "r");
	fscanf(fp, "%lf%d%d", &beta, &i, &period);
	for (i = 0; i < N; i++) {
		fscanf(fp, "%lf", &h[i]);
	}
	for (i = 0; i < N; i++) {
		fscanf(fp, "%d", &K[i]);
	}
	fscanf(fp, "%lf%d%d%d", &p, &UB, &LB, &D_len);
	for (i = 0; i < D_len; i++) {
		fscanf(fp, "%d", &D[i]);
	}
}

int main(int argc, const char *argv[])
{
	int X[3];
	init();
	DP_U1();
	DP_U2();
	while(scanf("%d%d%d", &X[0], &X[1], &X[2]) != EOF){
		printf("%d\t%d\t%d\t", X[0], X[1], X[2]);
		printf("%.2lf\t", MAX(get_U1(X), get_U2(X)));
	}
	return 0;
}
