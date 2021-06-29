#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils.hpp"

const std::map<std::string,std::string> Filenames({{"H1","H1(el, 0).dat"},
											{"He4","He4(el, 0).dat"},
											{"O16","O16(el, 0).dat"}});

const std::string result = "Solar(0).dat";

int main(){

	
	std::vector<double> Mgrid = grid(1000,2,100,1,0);
	
	std::vector<double> P(Mgrid.size());
	
	for(auto el: Filenames ){
		auto F = Function::LoadFunction1(el.second.c_str());
		std::cout << "element = " <<el << ": " << F.getXgrid().size();
		
		double N = ME.at(el.first);
		P  = P + apply_function<double,double>(Mgrid,[N,&F](double mk){
			return F(mk)*N*(1+mk)/(N+mk)/mk;
		});
	}
	
	std::ofstream out(result);
	out << "Mw\tP"<<std::endl;
	out << Function1<double>(Mgrid,P);
	
}