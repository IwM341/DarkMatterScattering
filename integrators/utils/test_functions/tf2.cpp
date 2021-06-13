#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>

int main(void){

	std::vector<double> X(10);
	std::vector<double> Y(10);
	
	for(size_t i=0;i<10;i++){
		X[i] = i;
		Y[i] = i;
	}
	auto f = [](double x,double y){return x*y;};
	Function2<double> F2(X,Y,f);
	std::cout << F2.toString() << std::endl;
	
	for(size_t i=0;i<120;i++){
		double x = 10*((double)std::rand() )/RAND_MAX;;
		double y = 10*((double)std::rand() )/RAND_MAX;;
		std::cout << "f(" << x << ", " << y << ") = " <<
				F2(x,y) <<" vs " << f(x,y) << std::endl;
	}	
	return 0;
}
