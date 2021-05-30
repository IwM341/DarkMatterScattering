#include <cstdio>
#include <iostream>
#include "integrator.hpp"
#include "phase4.hpp"
#include "cross_section3.hpp"

#include <fstream>
#include <string>

#include <time.h>

#define a1 1
#define b1 0
#define a2 1
#define b2 0

double delta(double x,double sigma){
	return exp(-x*x/(2*sigma*sigma))/sqrt(2*M_PI*sigma*sigma);
}
double M2phi(PhaseState3 out, PhaseState2 in){
	
	const double e2 = (4*M_PI/137.036);
	double p1q = out.P1*out.Q;
	double pq = in.P0*out.Q;
	
	
	double M2hi = 4*((a2*a2+b2*b2)*out.K1*in.K0 + (a2*a2-b2*b2)*in.K0.quad());
	double M2ee = ((a1*a1+b1*b1)*out.P1*in.P0 + (a1*a1-b1*b1)*in.P0.quad());
	
	double dM2ee = (a1*a1+b1*b1)*(in.P0-out.P1)*out.Q;
	

	double qcoeff = -(out.P1/p1q - in.P0/pq).quad();
	
	double part3 =  (a1*a1+b1*b1)*(pq-p1q)*(1/p1q-1/pq);
	
	return M2hi*(qcoeff*(M2ee+dM2ee)+part3)*e2;
}

template <typename Functor2>
void save_function2(Functor2 F,double x0,double x1,double y0,double y1,int Nx,int Ny,
					std::string filename = "save_function.dat",
					std::string xtitle = "x",std::string ytitle = "y"){
						
	double dx = (x1-x0)/(Nx-1);
	double dy = (y1-y0)/(Ny-1);
	
	std::ofstream out(filename.c_str(),std::ofstream::out);
	
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
	double mk = 0.1;
	double mp = 1;
	double k0 = 0.1;
	MatrixElementType M2 = M2phi;//[](PhaseState3,PhaseState2){return 1;};
	
	const double delta = 0.05;
	
	const double eps = delta*(E(mp,k0)+E(mk,k0))*(E(mk,k0)-mk)/(E(mk,k0)+E(mp,k0)-mk);
	
	
	auto F1 = dsigma_dk1_dcosTh_Q(mk,mp,k0,eps,M2,40,40);
	
	
	const double Ek0 = E(mk,k0);
	const double Ecm = Ek0 + E(mp,k0);
	
	const double k1p = ((1-eps/Ecm)*sqrt(k0*k0+eps*eps+2*eps*mk*mk/Ecm-2*eps*Ek0) 
							- eps*(Ek0-eps)/Ecm)/(1-2*eps/Ecm);
	
	const double k1m = ((1-eps/Ecm)*sqrt(k0*k0+eps*eps+2*eps*mk*mk/Ecm-2*eps*Ek0) 
							+ eps*(Ek0-eps)/Ecm)/(1-2*eps/Ecm);
	
	//std::cout << "km " << F1(k1m,0.5) << std::endl;
	//std::cout << "kp " << F1(k1p,0.5) << std::endl;
	
	save_function2(F1,0,k0,-1,1,1000,20,"testQ.dat","k'","cos(theta)");
	
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

