#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils.hpp"

#define SVAR(x) (std::string(#x) + std::string(" = ") + std::to_string(x))
#define PVAR(x) std::cout << std::string(#x) + std::string(" = ") + std::to_string(x) <<std::endl;

namespace std{
	string to_string(const string & S){
		return S;
	}
}

int main(void){
	double Vesc = 3.7336e-5;
	std::vector<double> Ugrid({0.01,0.03,0.05,0.063,0.0645,0.0646,0.0647,0.0648,0.0649,
										0.065,0.0710,0.08,0.1,0.3,0.7,1.0});
	
	double mk = 4;
	double mp = 16;
	
	auto M22 = MET0(mp+mk);
	auto M23 = MET0Q(mp+mk);
	
	std::cout << "eval SIn" <<std::endl;
	std::vector<double> SIn = apply_function<double,double>(Ugrid*U0,[Vesc,M22,mp,mk](double u){
							double v = sqrt(Vesc*Vesc+u*u);
							return sigmaC(mp,mk,v,Vesc,M22,10);
						});
						
	std::cout << "eval Srg" <<std::endl;
	
	std::vector<double> SRg = apply_function<double,double>(Ugrid*U0,[Vesc,M22,M23,mp,mk](double u){
							double v = sqrt(Vesc*Vesc+u*u);
							return sigmaC(mp,mk,v,Vesc,M22,M23,6,6,10,40,0.01);
						});
						
	std::cout << "eval Sel" <<std::endl;
	
	std::vector<double> Sel = apply_function<double,double>(Ugrid*U0,[Vesc,mp,mk](double u){
							double v = sqrt(Vesc*Vesc+u*u);
							return sigmaTfacor(mp,mk,v,Vesc,U0,0,CAPTURE,0,1);
						});
	
	
	
	std::vector<double> UgridND = Ugrid;
	
	std::ofstream out1("sin\\Sin.dat");
	out1 << "U\tSgm" << std::endl;
	out1 << Function1<double>(UgridND,SIn);
	
	std::ofstream out2("sin\\SRg.dat");
	out2 << "U\tSgm" << std::endl;
	out2 << Function1<double>(UgridND,SRg);
	
	std::ofstream out3("sin\\Sel.dat");
	out3 << "U\tSgm" << std::endl;
	out3 << Function1<double>(UgridND,Sel);
}