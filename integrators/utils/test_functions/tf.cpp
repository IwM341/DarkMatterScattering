#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>

int main(void){
	std::srand(time(NULL));

	size_t N = 100;
	std::vector<double> X(N);
	for(size_t i=0;i<N;i++)
		X[i] = ((double)std::rand() )/RAND_MAX;
	
	std::sort(X.begin(),X.end());
	
	double x = ((double)std::rand() )/RAND_MAX;
	size_t i = find_less(X,x);
	size_t i1 = i+1;
	if(i1 >= X.size())
		i1 = i;
	std::cout << "X.size() = " <<X.size() <<std::endl;
	std::cout << X <<std::endl;
	std::cout << X[i] << " < " << x << " < " <<X[i1] <<std::endl;
	
	
	Function1<double> F1(X,apply_function<double,double>(X,[](double x)->double{return x;}));
	auto F = F1.evald();
	
	/*
	std::ofstream out("filename.txt");
	for(int i=0;i<10;i++){
		x = ((double)std::rand() )/RAND_MAX;
		std::cout << "f(" << x <<") = "<<F(x) << " vs " << x << std::endl;   
	}
	*/	
	Function1<double> F2 = Function::LoadFunction1("filename.txt");
	std::cout << F2 << std::endl;
	return 0;
}
