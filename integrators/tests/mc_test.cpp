#include "mc.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>

double delta(double x,double sgm){
	return std::exp(-x*x/(2*sgm*sgm))/std::sqrt(2*M_PI*sgm*sgm);
}

double f(double x){
	return delta(x-0.5,0.01);
}
class my_gen:public AreaGenerator<double>{
	public:
	my_gen(){ srand(time(0));}
	virtual double RandomCoords() const{
		//double x1 = std::rand()*2.0/RAND_MAX-1;
		//double x2 = std::rand()*2.0/RAND_MAX-1;
		//return std::vector<double> ({x1,x2});
		return std::rand()*1.0/RAND_MAX;
	}
	virtual double AreaVolume() const{return 1;	}
};



int main(){
	
	
	MonteCarloIntegrator I;
	my_gen gen;
	std::cout<< I(f,&gen,100000) << std::endl;
	
	return 0;
	
}
