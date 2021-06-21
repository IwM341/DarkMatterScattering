#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>


#include "../utils.hpp"

int main(void){
	/*
	auto &EM = BodyModel::Instance();
	std::vector<double> Ugrid(100);
	std::vector<double> VescGrid(10);
	
	double Nu = Ugrid.size();
	double Nve = VescGrid.size();
	
	for(size_t i=0;i<Ugrid.size();i++)
		Ugrid[i] = 4*i*U0/Nu;
	for(size_t i=0;i<VescGrid.size();i++)
		VescGrid[i] = EM.VeMin() + i*(EM.VeMax() - EM.VeMin())/Nve;
	
	auto SI = SigmaInelastic(2.0,"Fe",Ugrid,VescGrid);
	
	std::cout << SI <<std::endl;
	
	auto Su = IntegrateSigma(SI,"Fe");
	

	std::ofstream out(std::string("SigmaU_m(") + std::to_string(2.0) + ").dat");
	
	out << "U\tSigma\n" << Su <<std::endl;
	
	std::cout <<"*/
	
	double cpt = sigmaTfacor(10,1,1,0.1,1,0.0,CAPTURE,2,1000000);
	double esc =  sigmaTfacor(10,1,0.5,1,1,0.1,ESCAPE,2,1000000);
	
	std::cout << cpt << std::endl;
	std::cout << esc << std::endl;
	return 0;
}
