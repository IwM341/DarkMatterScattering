#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils.hpp"

#define Tp "el"
#define tp 0

const std::map<std::string,std::string> Filenames({{"Fe",std::string("Fe(") + Tp + ", "+ std::to_string(tp) + ").dat"},
											{"Ni",std::string("Ni(") + Tp + ", "+ std::to_string(tp) + ").dat"},
											{"Mg",std::string("Mg(") + Tp + ", "+ std::to_string(tp) + ").dat"},
											{"Si",std::string("Si(") + Tp + ", "+ std::to_string(tp) + ").dat"},
											{"O",std::string("O(") + Tp + ", "+ std::to_string(tp) + ").dat"}});

const std::string result = std::string(Tp) + "(" + std::to_string(tp)+ ").dat";

int main(){

	
	std::vector<double> Mgrid = grid(1000,2,100,1,0);
	
	std::vector<double> P(Mgrid.size());
	
	for(auto el: Filenames ){
		auto F = Function::LoadFunction1(el.second.c_str());
		std::cout << F.getXgrid().size();
		double N = ME.at(el.first);
		P  = P + apply_function<double,double>(Mgrid,[N,&F](double mk){
			return F(mk)*N*(1+mk)*(1+mk)/(N+mk)/(N+mk)/mk;
		});
	}
	
	std::ofstream out(result);
	out << "Mw\tP"<<std::endl;
	out << Function1<double>(Mgrid,P);
	
}