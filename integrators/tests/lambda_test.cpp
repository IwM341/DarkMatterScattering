#include <iostream>
#include <functional>

int main(void){
	
	auto F = [](double x){return x;};
	std::cout << F(1) << std::endl;
	
	return 0;
}
