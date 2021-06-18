#include <cstdio>
#include <iostream>
#include <functional>
int main(void){
	std::function<double(double)> X = [](double x){
		int N = 10000;
		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){
				x = x*(1.0 + x*(i+j)/(2.0*N));
			}
		}
		return x;
	};
	
	//#pragma omp parallel for 
	for(int i=0;i<64;i++)
		std::cout << X(i/64.0) << "\t";
	std::cout << std::endl;
}
