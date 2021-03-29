#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include "phase4.hpp"

double F1(double (*f)(double)){
	return f(8);
}

int main(){
	vec3 x(1,2,3);
	std::cout<< (std::string)(2*(vec4(vec3(0,0,1),4) - vec4(vec3(0,0,0),4))) << std::endl;
	return 0;
}
