#include <cstdio>
#include <iostream>
#include "../integrators/lambda_int.hpp"
#include "../phase4.hpp"
#include "../sections/cross_section3.hpp"

#include <fstream>
#include <string>

#include <time.h>
#define a1 1
#define b1 0
#define a2 1
#define b2 0

double M2phi(PhaseState3 out, PhaseState2 in){
	
	double p1q = out.P1*out.Q;
	double pq = in.P0*out.Q;
	
	
	double M2hi = 4*((a2*a2+b2*b2)*out.K1*in.K0 + (a2*a2-b2*b2)*in.K0.quad());
	double M2ee = ((a1*a1+b1*b1)*out.P1*in.P0 + (a1*a1-b1*b1)*in.P0.quad());
	
	double dM2ee = (a1*a1+b1*b1)*(in.P0-out.P1)*out.Q;
	
	double qcoeff = -(out.P1/p1q - in.P0/pq).quad();
	
	double part3 =  (a1*a1+b1*b1)*2*(pq-p1q)*(1/p1q-1/pq);
	
	return M2hi*(qcoeff*(M2ee+dM2ee)+part3);
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

int main(void){
	double mk = 1;
	double mp = 0.1;
	double k0 = 0.001;
	MatrixElementType M2 = M2phi;//[](PhaseState3,PhaseState2){return 1;};
	auto F1 = dsigma_dk1_dcosTh(mk,mp,k0,M2,80,80);
	
	save_function2(F1,0,k0,-1,1,40,40,"sigmas1.dat","k'","cos(theta)");
}

