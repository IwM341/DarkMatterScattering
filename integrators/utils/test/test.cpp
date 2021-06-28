#include <cstdio>
#include <iostream>
#include "integrator.hpp"
#include "phase4.hpp"
#include "cross_section3.hpp"
#include "matrix_element.hpp"
#include <utility> 
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>

#define a1 1
#define b1 0


template <typename Functor2>
void save_function2(Functor2 F,double x0,double x1,double y0,double y1,int Nx,int Ny,
					std::string filename = "save_function.dat",
					std::string xtitle = "x",std::string ytitle = "y"){
						
	double dx = (x1-x0)/(Nx-1);
	double dy = (y1-y0)/(Ny-1);
	
	std::ofstream out(filename.c_str(),std::ofstream::out);
	//out <<std::setprecision(17);
	out << ytitle+"\\"+xtitle;
	for(int j =0;j<Nx;j++){
		out << '\t' << x0 + j*dx;
	}
	out << '\n';
	for(int i =0;i<Ny;i++){
		out << y0 + i*dy;
		for(int j=0;j<Nx;j++){
			out << '\t' << F(x0 + j*dx,y0 + i*dy);
		}
		out << '\n';
	}
}

const double pbarn_to_GeV = 2.56818998849288e-09;
const double GeV_to_pbarn = 1/2.56818998849288e-09;

int main(void){
	double mk = 0.5;
	double mp = 1;
	double k0 = 0.001;
	double Ep = E(mp,k0);
	double Ek = E(mk,k0);
	
	MatrixElementType23 M23 = MET1Q(Ep+Ek,k0);//[](PhaseState3,PhaseState2){return 1;};
	MatrixElementType22 M22 = MET1(Ep+Ek,k0);
	
	
	
	auto F1 = dsigma_d3k1(mk,mp,k0,M23,40,40);
	auto F2 = dsigma_d3k1_NR(mk,mp,k0,M22);
	
	save_function2(F1,0,k0,-1,1,20,20,"ordinary.dat","k'","cos(theta)");
	save_function2(F2,0,k0,-1,1,20,20,"simplify.dat","k'","cos(theta)");
	
	/*
	double cosTh = 0.5;
	//double k1 = 0.5*k0;
	double Ecm = E(mk,k0)+E(mp,k0);
	//double E1 = E(mk,k1);
	double E0 = E(mk,k0);
	
	
	std::cout << integrateAB([E0,k0,Ecm,cosTh,F1](double k1){
			return F1(k1,cosTh)*2*GeV_to_pbarn /*4*k0*Ecm;
		},0,k0,100) << std::endl;
	
	std::cout << integrateAB([E0,k0,Ecm,cosTh,F1,mk](double k1){
			double E1 = E(mk,k1);
			return k1*k1*Ecm*(E0-E1)/(2*E1*phase_2pi3*((Ecm-E1)*(Ecm-E1) - k1*k1));
		},0,k0,100) << std::endl;
	/*
	//std::cout <<"dsigma/domega = " << F1(k1,cosTh)*2*4*k0*Ecm << std::endl;
	//std::cout << k1*k1*Ecm*(E0-E1)/(2*E1*phase_2pi3*((Ecm-E1)*(Ecm-E1) - k1*k1)) << std::endl;
	//std::cout << F1(k1,cosTh)*2*4*k0*Ecm -  k1*k1*Ecm*(E0-E1)/(2*E1*phase_2pi3*((Ecm-E1)*(Ecm-E1) - k1*k1))  << std::endl;
	/**/
	
	
	return 0;
}

