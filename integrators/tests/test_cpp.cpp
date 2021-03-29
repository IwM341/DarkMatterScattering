#include <cstdlib>
#include <iostream>

double f(double x,double y){
	return x+y;
}

#define INT(F,a,b) (F(a)+F(b))/2

double g(double x){
	#define F(y) f(x,y)
	return INT(F,0,1);
}

int main(void){
	int i=0;
	{
		int i = 1;
		std::cout << i <<std::endl;
	}
	std::cout << i <<std::endl;
	std::cout << g(1) << std::endl;
	return 0;
}
