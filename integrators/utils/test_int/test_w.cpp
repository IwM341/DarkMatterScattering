#include <cstdio>
#include <iostream>
#include "integrator.hpp"
#include "cross_section3.hpp"
#include <cmath>

#define ABS(X) ((X) > 0 ? X : -X)
double compare(double x,double y){
	return log2(ABS(2*(x-y)/(x+y)));
}
int main(void){
	
	double x1 = 0.001;
	double y1 = 6.6666706666695238117e-7;
	
	double x2 = 0.0001;
	double y2 = 6.6666667066666669524e-9;
	
	double x3 = 0.01;
	double y3 = 0.66670666952403176422e-4;
	
	double x4 = 0.1;
	double y4 = 0.67069546215116127145e-2;
	
	double x5 = 0.3;
	double y5 = 0.63464028020744769827e-1;
	
	double x6 = 0.5;
	double y6 = 0.19722457733621938279;
	
	double x7 = 0.7;
	double y7 = 0.47800150769729484122;
	
	std::cout <<"comparation:"<< std::endl;
	std::cout << compare(W(x1),y1) <<std::endl;
	std::cout << compare(W(x2),y2) <<std::endl;
	std::cout << compare(W(x3),y3) <<std::endl;
	std::cout << compare(W(x4),y4) <<std::endl;
	std::cout << compare(W(x5),y5) <<std::endl;
	std::cout << compare(W(x6),y6) <<std::endl;
	std::cout << compare(W(x7),y7) <<std::endl;
	return 0;
}
