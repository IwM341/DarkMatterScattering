#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils_S.hpp"

const std::map<std::string,std::string> Filenames({{"H1","H1(1).dat"},
											{"He4","He4(1).dat"},
											{"O16","O16(1).dat"}});

const std::string result = "Solar(1).dat";

int main(){

	
	std::vector<double> Mgrid = grid(1000,2,100,1,0);
	
	std::vector<double> P(Mgrid.size());
	
	for(auto el: Filenames ){
		auto F = Function::LoadFunction1(el.second.c_str());
		std::cout << "element = " <<el.first << ": " << F.getXgrid().size();
		
		double N = ME.at(el.first);
		P  = P + apply_function<double,double>(Mgrid,[N,&F](double mk){
			return F(mk)*N*(1+mk)*(1+mk)/(N+mk)/(N+mk)/mk;
		});
	}
	
	std::ofstream out(result);
	out << "Mw\tP"<<std::endl;
	out << Function1<double>(Mgrid,P);
	
}