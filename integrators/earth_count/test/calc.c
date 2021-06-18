#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "../utils.hpp"

int main(void){
	auto &EM = BodyModel::Instance();
	Function2<double> Sigma(std::vector<double>({0.,1.,2.,3.,4.,5.}),
							EM["Vesc"],[](double x,double y){return x*y;});
	
	std::cout << Sigma <<std::endl;
	
	auto Su = IntegrateSigma(Sigma,"H1");
	
	
	
	std::cout << Su <<std::endl;
	
	return 0;
}
