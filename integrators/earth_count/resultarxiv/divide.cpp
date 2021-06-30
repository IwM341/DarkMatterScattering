#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils.hpp"

int main(int argc,char**argv){
	if(argc < 4){
		std::cout << "exprct [inelastic] [elastic] [out]" <<std::endl;
		return 0;
	}
	auto IN = Function::LoadFunction1(argv[1]);
	auto EL = Function::LoadFunction1(argv[2]);
	
	Function1<double> DIV = IN;
	DIV.insert(EL);
	
	
	DIV = Function1<double>(DIV.getXgrid(),[&IN,&EL](double mass){return EL(mass);});
	
	//std::cout << DIV << std::endl;
	
	std::ofstream out(argv[3]);
	out << "mass\tdivision"<<std::endl;
	out << DIV;
	
	
	return 0;
}