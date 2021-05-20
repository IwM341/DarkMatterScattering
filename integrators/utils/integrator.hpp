#ifndef INT_LAMBDA
#define INT_LAMBDA

#include "omp.h"
#include <functional>
#include <utility>
#include <cmath>


						
						
						
extern inline double integrateAB(std::function<double(double)> f,double a,double b,int N){
	
	double sum = 0.0;
	double sum5;
	double xi = 0;
	double h = (b-a)/N;
	#pragma omp parallel for reduction(+:sum) private(sum5,xi) shared(N,h)
	for (int i = 0; i < N; i++) {
		xi = a + i*h;
		sum5 = f(xi +0.5*h);/*
		for(int j = 0;j<5;j++)
			sum5 += Wh[j]*f(xi + h*Xh[j]);*/
		sum += sum5*h;
	}
	return sum;
}

const double Xh[5] = {0.5, 0.5+0.5*.53846931010568309105, 0.5-0.5*.53846931010568309105,
					0.5+0.5*.90617984593866399282, 0.5-0.5*.90617984593866399282};
const double Wh[5] = {.28444444444444444444, .23931433524968323402, .23931433524968323402,
					.11846344252809454376, .11846344252809454376};

extern inline double integrateAB5(std::function<double(double)> f,double a,double b,int N){

	double sum = 0.0;
	double sum5;
	double xi = 0;
	double h = (b-a)/N;
	#pragma omp parallel for reduction(+:sum) private(sum5,xi) shared(N,h)
	for (int i = 0; i < N; i++) {
		sum5 = 0.0;
		xi = a + i*h;
		for(int j = 0;j<5;j++)
			sum5 += Wh[j]*f(xi + h*Xh[j]);
		sum += sum5*h;
	}
	return sum;
}
extern inline double integrateLogAB5(std::function<double(double)> f,double a,double b,double eps,int N){
	const double L = b-a;
	
	const double T = pow(eps/(L+eps),1/N);
	
	double sum = 0.0;
	double sum5;
	
	
	double h = (b-a)/N;
	#pragma omp parallel for reduction(+:sum) private(sum5,xi) shared(N,h)
	for (int i = 0; i < N; i++) {
		sum5 = 0.0;
		const double xi = a - eps + (L+eps)*pow(T,i+1);
		const double xi1 = a - eps + (L+eps)*pow(T,i);
		const double h = xi1 - xi;
		for(int j = 0;j<5;j++)
			sum5 += Wh[j]*f(xi + h*Xh[j]);
		sum += sum5*h;
	}
	return sum;
}				
#endif
