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
	
	
	
	Function1<double> F1(X,apply_function<double,double>(X,[](double x)->double{return x;}));
	
	/*
	std::ofstream out("filename.txt");
	for(int i=0;i<10;i++){
		x = ((double)std::rand() )/RAND_MAX;
		std::cout << "f(" << x <<") = "<<F(x) << " vs " << x << std::endl;   
	}
	*/	
	auto f = (std::function<double(double)>)F1;
	Function1<double> F2 = Function::LoadFunction1("filename.txt");
	std::cout << f(3.3) << std::endl;
	return 0;
}
