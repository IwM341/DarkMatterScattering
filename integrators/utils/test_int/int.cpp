#include <cstdio>
#include <iostream>
#include "integrator.hpp"
#include <cmath>
int main(void){
	double eps = 0.001;
	int N = 10;
	std::cout << integrateRightLogAB5([eps](double x){return 1.0/(1.+eps-x);},
	0.,1.,eps,N) <<" vs "<< integrateAB5([eps](double x){return 1.0/(1.+eps-x);},
	0.,1.,N) << " and " <<
	log((1.0+eps)/eps) << std::endl;
	return 0;
}
