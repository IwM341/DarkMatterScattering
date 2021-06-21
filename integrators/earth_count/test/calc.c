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
	
	std::cout << sigmaTfacor(10,10,1,1,1,1,ESCAPE,0,1000000) << std::endl;
	return 0;
}
