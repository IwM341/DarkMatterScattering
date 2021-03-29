//#include "int_macros.h"
#include "lambda_int.hpp"

#include <time.h>
#include <iostream>


typedef std::function<double(double)> ftype; 

#define M 3000

double f(double x,double y){
	return x+y;
}



double intL(){
	return integrateAB((ftype)[](double x){
		return integrateAB((ftype)[x](double y){return f(x,y);},M,0,1);
	},M,0,1);
}

int main(){
	
	double t = clock();
	
	std::cout << intL() << ", time = " << (clock()-t) << std::endl;
	
	return 0;
}
