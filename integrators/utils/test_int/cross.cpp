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
	double mk = 0.7;

	double v_ls = 0.001;
	
	double v_esc = v_ls/3;
	

	

	//std::cout << sigmaC(mp,mk,v_ls,v_esc,M22,M23,6,6,6,100,0.01) <<std::endl;
	
	const size_t M = 100;
	
	std::vector<double> mks(M);
	std::vector<double> sgms(M);
	auto M22 = MET1(mp+mk,v_ls);
	auto M23 = MET1Q(mp+mk,v_ls);
	
	//auto M22 = ScalarElement(mp,mk);
	//auto M23 = ScalarElementQ(mp,mk);
	
	std::cout << sigmaC(mp,mk,v_ls,v_esc,M22,M23,4,4,4,40,0.001) <<std::endl;
	std::cout << sigmaC(mp,mk,v_ls,v_esc,M22,M23,4,4,4,40,0.01) <<std::endl;
	std::cout << sigmaC(mp,mk,v_ls,v_esc,M22,M23,4,4,4,40,0.1) <<std::endl;
	
	/*
	#pragma omp parallel for shared(mks,sgms)
	for(size_t i = 0;i<M;i++){
		auto M22 = MET0(mp+mk);
		auto M23 = MET0Q(mp,mk);
		double mk = 0.4+0.01*i;
		mks[i] = mk;
		sgms[i] = (sigmaC(mp,mk,v_ls,v_esc,M22,M23,2,2,2,50,0.01));
	}
	
	std::cout << mks <<std::endl;
	std::cout << sgms <<std::endl;
	*/
	return 0;
}
