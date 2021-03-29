#ifndef INT_MACROS
#define INT_MACROS

#include "omp.h"

#define X1 0.5
#define X2 0.76923465505284154553
#define X3 0.23076534494715845447
#define X4 0.95308992296933199641
#define X5 0.04691007703066800359

#define W1 0.28444444444444444444
#define W2 0.23931433524968323402
#define W3 0.23931433524968323402
#define W4 0.11846344252809454376
#define W5 0.11846344252809454376

/*
#define INTEGRATE_OMP(F,a,b,N,RES) {\
	double __sum##RES = 0.0;\
	double __xi##RES = 0;\
	double __h##RES = (b-a)/N;\
	#pragma omp parallel for reduction(+:__sum##RES) private(__xi##RES)\
	for (int i = 0; i < N; i++) {\
		_xi##RES = a + i*h;\
		__sum##RES += h*(W1*F(xi + h*X1) +W2*F(xi + h*X2) + W3*F(xi + h*X3) + \
							W4*F(xi + h*X4) + W5*F(xi + h*X5));\
	}\
	RES = __sum##RES;\
}


#define INTEGRATE(F,a,b,N,RES) {\
	double __sum##RES = 0.0;\
	double __xi##RES = 0;\
	double __h##RES = (b-a)/N;\
	for (int i = 0; i < N; i++) {\
		_xi##RES = a + i*h;\
		__sum##RES += h*(W1*F(xi + h*X1) +W2*F(xi + h*X2) + W3*F(xi + h*X3) + \
							W4*F(xi + h*X4) + W5*F(xi + h*X5));\
	}\
	RES = __sum##RES;\
}
*/
extern inline double integrateOMP(double(*f)(double),double a,double b,double N){
	double sum = 0.0;
	double xi = 0;
	double h = (b-a)/N;
	#pragma omp parallel for reduction(+:sum) private(ix)
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum += (W1*f(xi + h*X1) +W2*f(xi + h*X2) + W3*f(xi + h*X3) + 
							W4*f(xi + h*X4) + W5*f(xi + h*X5));
	}
	return h*sum;
}

extern inline double integrate(double(*f)(double),double a,double b,double N){
	double sum = 0.0;
	double xi = 0;
	double h = (b-a)/N;
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum += (W1*f(xi + h*X1) +W2*f(xi + h*X2) + W3*f(xi + h*X3) + 
							W4*f(xi + h*X4) + W5*f(xi + h*X5));
	}
	return h*sum;
}
#endif
