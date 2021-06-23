#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>


#include "../utils.hpp"
const std::vector<std::string> EarthElements({"Fe","O","Si","Mg","Ni","Ca","Al"}); 
/*
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
	
	
	double cpt = sigmaTfacor(10,1,1,0.1,1,0.0,CAPTURE,2,1000000);
	double esc =  sigmaTfacor(10,1,0.5,1,1,0.1,ESCAPE,2,1000000);
	
	std::cout << cpt << std::endl;
	std::cout << esc << std::endl;
	return 0;
}
*/

std::vector<double> grid(size_t N,double Umin,double Umax,double q,double eps){
	std::vector<double> __grid__(N);
	
	double et = pow(eps,1/q);
	double Z = pow(1+eps,1/q)-et; 
	
	for(size_t i=0;i<N;i++){
		double t = (double)i/N;
		__grid__[i] = (pow(t*Z+et,q)-eps)*(Umax-Umin) + Umin;
	}
	return __grid__;
}

void save_F2(const Function1<double> &F2,
				const std::string path,const std::string element,
				int isElastic,double mk,int sigmaType){
	std::string elastictype = "In";
	if(isElastic)
		elastictype = "El";
	std::ofstream out(path + elastictype + "(" + element +
		", " + std::to_string(mk) + ", " +  std::to_string(sigmaType) + ").dat");
	out << "U\tSigma\n" << F2 <<std::endl;
}

int main(void){
	auto &EM = BodyModel::Instance();
	std::vector<double> UgridInelastic = grid(40,0.0,U0*3,2,0.01);
	std::vector<double> UgridElastic = grid(40,0.0,U0/2,4,0.01);
	std::vector<double> VescGrid = grid(16,EM.VeMin(),EM.VeMax(),1,0.0);
	
	std::vector<double> massGrid = grid(8,2,5,1,0.0);
	
	//std::cout << Ugrid <<std::endl;
	//std::cout << VescGrid <<std::endl;
	
	double Xi = 0.5;
	for(int sigma_type =0;sigma_type<3;sigma_type++){
		for(auto el:EarthElements){
			std::vector<double> CIn;
			std::vector<double> CEl;
			for(double mass:massGrid){
				Function2<double> SI = SigmaInelastic(mass,el,UgridInelastic,VescGrid,sigma_type);
				Function2<double> SE = SigmaElastic(mass,el,UgridElastic,VescGrid,0,sigma_type,10000);
				
				Function1<double> SuI = IntegrateSigma(SI,el);
				Function1<double> SuE = IntegrateSigma(SE,el);
				
				save_F2(SuI,"results\\",el,0,mass,sigma_type);
				save_F2(SuE,"results\\",el,1,mass,sigma_type);
				
				CIn.push_back(C_ND(SuI,Xi*U0));
				CEl.push_back(C_ND(SuE,Xi*U0));
			}
			Function1<double> FIn(massGrid,CIn,true);
			Function1<double> FEl(massGrid,CEl,true);
			
			std::ofstream outEl(std::string("Cel(")+std::to_string(sigma_type) + 
										", " + el + ").dat");
			std::ofstream outIn(std::string("Cin(")+std::to_string(sigma_type) + 
										", " + el + ").dat");
										
			outEl << "mass\tC\n" << FEl <<std::endl;
			outIn << "mass\tC\n" << FIn <<std::endl;
		}
	}
	return 0;
}



/*
int main(void){
	auto &EM = BodyModel::Instance();
	for(auto el:EarthElements){
		std::cout << el<< std::endl;
		std::cout << EM[el].size() << std::endl;
		std::cout << ME.at(el) << std::endl<< std::endl;
	}

	
	return 0;
}
*/
