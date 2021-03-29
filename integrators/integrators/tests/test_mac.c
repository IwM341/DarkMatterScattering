
#include "int_macros.h"
#include <stdio.h>
#include <time.h>


#define M 3000

double f(double x,double y){
	return x+y;
}

double intM(){

	inline double Sc(double x)
	{
		inline double Fc(double y){
			return f(x,y);
		}
		return integrate(Fc,0,1,M);
	}
	return integrate(Sc,0,1,M);
}

int main(void){
	double  t =  clock();
	double x = intM();
	printf("%lf, time = %lf",x,clock()-t);
	return 0;
}
