#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <chrono>
//#include <time.h>

auto Functor(double N){
	return [N](){
			double sum = 0;
			for(int i=0;i<N;i++){
				sum += i;
			}
			return sum;
		};
}

int main(void){
	const int M = 80;
	double a[M];
	int N = 10000000;
	
	double now1 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	
	std::cout << "now = " << now1 <<std::endl;
	
	#pragma omp parallel for 
	for(int i=0;i<M;i++){
		a[i] = Functor(N+i)();
	}
	
	double after = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	
	std::cout << "after = " << after <<std::endl;
	
	double delta = after - now1;
	
	std::cout << "delta = " << delta <<std::endl;
	
	
	for(int i=0;i<7;i++){
		std::cout << a[i] << ", ";
	}
	std::cout << a[7] << std::endl;
	
	return 0;
}
