#include <cstdio>
#include <iostream>
#include <ostream>
#include "integrator.hpp"
#include "cross_section3.hpp"
#include "matrix_element.hpp"
#include <vector>
//#include <omp.h>
#include <cmath>

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& V)
{
	os << "vector[";
	for(size_t i = 0;i<V.size()-1;i++){
		os << V[i] << ", ";
	}
	if(V.size() != 0)
		os << V[0];
	os << "]";
    return os;
}


int main(void){
	
	double mp = 1;
	double mk = 0.1;

	double v_ls = 0.001;
	
	double v_esc = v_ls/3;
	

	

	//std::cout << sigmaC(mp,mk,v_ls,v_esc,M22,M23,6,6,6,100,0.01) <<std::endl;
	
	const double M = 20;
	
	std::vector<double> mks(M);
	std::vector<double> sgms(M);
	
	#pragma omp paralel for
	for(size_t i = 0;i<M;++i){
		auto M22 = ScalarElement(mp,mk);
		auto M23 = ScalarElementQ(mp,mk);
		double mk = 0.4+0.01*i;
		mks[i] = mk;
		sgms[i] = (sigmaC(mp,mk,v_ls,v_esc,M22,M23,6,6,6,100,0.01));
	}
	
	std::cout << mks <<std::endl;
	std::cout << sgms <<std::endl;
	
	return 0;
}
