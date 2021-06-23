#ifndef CROSSSEC_H
#define CROSSSEC_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include "integrator.hpp" 
#include "omp.h" 
#include "phase4.hpp"
#include <functional>
#include "cstdlib"

#define  phase_2pi1   6.2831853071795864770
#define  phase_2pi2   39.478417604357434476
#define  phase_2pi3   248.05021344239856142
#define  phase_2pi4   1558.5454565440389959
#define  phase_2pi5   9792.6299131290065050

#ifdef FAST
#define INTEGRATOR integrateAB5
#else
#define INTEGRATOR integrateAB
#endif

extern inline double fase_volume_density(double k1,double Ecm,double Ek0,double Ek1,double cosTh1){
	double q_part = Ecm - Ek1 - k1*cosTh1;
	return ( k1* k1*Ecm*(Ek0-Ek1))/(8*q_part*q_part*phase_2pi4*Ek1);
}

//pase volume parameters: dk' d cosTheta d cosTheta1 d phi1

//integration through: d cosTheta1 d phi1


class PhaseState3{
public:
	vec4 K,P,Q;
	PhaseState3(const vec4 &K,const vec4 & P,const vec4 & Q):K(K),P(P),Q(Q){}
	PhaseState3(){}
	

};
class PhaseState2{
public:
	vec4 K,P;
	PhaseState2(const vec4 & K,const vec4 & P):K(K),P(P){}
	PhaseState2(){}
	

};



std::ostream& operator<<(std::ostream& os, const PhaseState3& V3){
	os << "k1 = " << V3.K.toString() << ", p1 = "<<V3.P.toString()<<",q = " <<V3.Q.toString();
	return os;
}

std::ostream& operator<<(std::ostream& os, const PhaseState2& V2){
	os << "k0 = " << V2.K.toString() << ", p0 = "<<V2.P.toString();
	return os;
}

typedef std::function<double(const PhaseState3 &out,const PhaseState2 &in)> MatrixElementType23;
typedef std::function<double(const PhaseState2 &out,const PhaseState2 &in)> MatrixElementType22;

typedef std::function<double(double k1,double cosTh)> dsigma_dk1_dcosTh_type;
typedef std::function<double(double k1,double cosTh)> dsigma_d3k1_type;

extern inline dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh(double mk,double mp,double k0,MatrixElementType23 M2,
								double Nphi ,double  NTh1){ 
	double Ek0 = E(mk,k0);
	double Ep0 = E(mp,k0);
	
	double Ecm = Ek0+Ep0;
	
	
	PhaseState2 St2(vec4(Ek0,0,0,k0),vec4(Ep0,0,0,-k0));
	
	return [mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,Nphi,NTh1](double k1,double cosTh){
		return integrateAB5([mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,k1,cosTh,Nphi,NTh1](double Theta1){
				return integrateAB5([mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,k1,cosTh,Theta1,Nphi,NTh1](double phi1){
						
						
						double cosTh1 = cos(Theta1);
						double sinTh1 = sin(Theta1);
						
						double Ek1 = E(mk,k1);
						
						double deltaE = (k0-k1)*(k0+k1)/(Ek0+Ek1);

						PhaseState3 St3;
						{
							double sinTh = sqrt(1-cosTh*cosTh);
							
							vec3 e3(sinTh,0,cosTh);
							vec3 e1(cosTh,0,-sinTh);
							vec3 e2(0,1,0);
							
							vec3 vk1 = k1*e3;
							
							
							double q =  (Ecm*(deltaE))/(Ecm- Ek1 + k1*cosTh1);
							
							vec3 vq = q*(cosTh1*e3+sinTh1*sin(phi1)*e2+sinTh1*cos(phi1)*e1);
							
							
							St3.K = vec4(Ek1,vk1);
							St3.Q = vec4(q,vq);
							St3.P = vec4(-vk1-vq,mp);
							
							/*
							#define COS(a,b) (a.vecPart()*b.vecPart()/a.vecPart().norm()/b.vecPart().norm())
							if(phi1<2 && phi1 > 1 && Theta1>1 && Theta1 < 2 && cosTh>0.4 & cosTh <0.6 && k1 > 0.5*k0){
								std::cout << "q = " << q << ", vq = " << vq.norm() << std::endl;
								std::cout << "phi ="<<phi1<<",theta1 = "<< Theta1 <<", cosTh = "<<cosTh<<", k1 = "<< k1 <<std::endl;
								std::cout << cosTh << "vs" << COS(St3.K,St2.K) << std::endl;
								std::cout << cosTh1 << "vs" << COS(St3.K,St3.Q) << std::endl;
								std::cout << "St2 = "<<St2<<std::endl;
								std::cout << "St3 = "<<St3<<std::endl;
								std::cout << "Ecm = " << Ecm << std::endl;
								std::cout << "k1^2,p1^2,q^2 = " << 
												St3.K.quad() << ", " <<
												St3.P.quad() << ", " << 
												St3.Q.quad() << std::endl;
								std::cout <<(St3.K+St3.Q+St3.P) <<std::endl;
								std::cout <<log ((St3.K+St3.Q+St3.P-St2.K - St2.P).t) <<std::endl;
								std::cout <<"Ecm - K'-P'-Q"  << Ecm - (St3.K+St3.Q+St3.P).t << std::endl;
								std::cout << std::endl;
							}/**/
						}
						
						double q_part = Ecm - Ek1 + k1*cosTh1;
						
						
						
						//std::cout<<deltaE<<" vs " << Ek0-Ek1 <<std::endl;
						
						double phase_density = 
								( k1*k1*Ecm*(deltaE))*sinTh1/(8*q_part*q_part*phase_2pi4*Ek1);
						double factor_p1p2 = 1/(4*k0*Ecm);
						
						if(q_part <0)
							std::cout<< Ecm << ", " << Ek1 << ", k1 =" << k1 << std::endl;
						
						return M2(St3,St2)*phase_density*factor_p1p2;
						
					},0,2*M_PI,Nphi);
			},0,M_PI,NTh1);
		};
	
}

extern inline dsigma_d3k1_type dsigma_d3k1(double mk,double mp,double k0,MatrixElementType23 M2,
								double Nphi ,double  NTh1){ 
	double Ek0 = E(mk,k0);
	double Ep0 = E(mp,k0);
	
	double Ecm = Ek0+Ep0;
	
	
	PhaseState2 St2(vec4(Ek0,0,0,k0),vec4(Ep0,0,0,-k0));
	double factor_p1p2 = 1/(4*k0*Ecm);
	return [mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,Nphi,NTh1,factor_p1p2](double k1,double cosTh){
		double sinTh = sqrt(1-cosTh*cosTh);
		vec3 e3(sinTh,0,cosTh);
		vec3 e1(cosTh,0,-sinTh);
		vec3 e2(0,1,0);
		
		vec3 vk1 = k1*e3;
		double Ek1 = E(mk,k1);
		double deltaE = (k0-k1)*(k0+k1)/(Ek0+Ek1);
		
		return integrateAB5([mk,mp,k1,Ek1,deltaE,factor_p1p2,vk1,e1,e2,e3,M2,Ecm,St2,Nphi](double Theta1){
				double cosTh1 = cos(Theta1);
				double sinTh1 = sin(Theta1);
				double q_part = Ecm - Ek1 + k1*cosTh1;
				double q =  (Ecm*(deltaE))/(Ecm- Ek1 + k1*cosTh1);
				double phase_density = 
								Ecm*deltaE*sinTh1/(8*q_part*q_part*phase_2pi5*Ek1);
								
				
				return integrateAB5([q,cosTh1,sinTh1,factor_p1p2,e1,e2,e3,vk1,St2,mp,Ek1,M2,phase_density](double phi1){				
						
						vec3 vq = q*(cosTh1*e3+sinTh1*sin(phi1)*e2+sinTh1*cos(phi1)*e1);
						PhaseState3 St3(vec4(Ek1,vk1),vec4(-vk1-vq,mp),vec4(q,vq));
						
						return M2(St3,St2)*phase_density*factor_p1p2;
						
					},0,2*M_PI,Nphi);
			},0,M_PI,NTh1);
		};
	
}

#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )

extern inline dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh_Q(double mk,double mp,double k0,double q_min,MatrixElementType23 M2,
								double Nphi ,double  NTh1){ 
	const double Ek0 = E(mk,k0);
	const double Ep0 = E(mp,k0);
	
	const double Ecm = Ek0+Ep0;
	
	
	const PhaseState2 St2(vec4(Ek0,0,0,k0),vec4(Ep0,0,0,-k0));
	
	const double k1p = ((1-q_min/Ecm)*sqrt(k0*k0+q_min*q_min+2*q_min*mk*mk/Ecm-2*q_min*Ek0) 
							- q_min*(Ek0-q_min)/Ecm)/(1-2*q_min/Ecm);
	
	const double k1m = ((1-q_min/Ecm)*sqrt(k0*k0+q_min*q_min+2*q_min*mk*mk/Ecm-2*q_min*Ek0) 
							+ q_min*(Ek0-q_min)/Ecm)/(1-2*q_min/Ecm);
	
	//std::cout<<"k1m = " << k1m << ",k1p = " << k1p << std::endl;
	//std::cout<<"cos_m = " << ( (Ecm-q_min)*(k0*k0-k1m*k1m)/(Ek0+E(mk,k1m)) - q_min*Ep0)/(q_min*k1m) << std::endl;
	//std::cout<<"cos_p = " << ( (Ecm-q_min)*(k0*k0-k1p*k1p)/(Ek0+E(mk,k1p)) - q_min*Ep0)/(q_min*k1p) << std::endl;

	return [q_min,mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,Nphi,NTh1,k1p,k1m](double k1,double cosTh){
		
		const double Ek1 = E(mk,k1);
		const double cosTh1_max = ( (Ecm-q_min)*(k0*k0-k1*k1)/(Ek0+Ek1) - q_min*Ep0)/(q_min*k1);
		
		const double sinTh = sqrt(1-cosTh*cosTh);
		
		const double Theta1_max = acos(MAX(-1,MIN(1,cosTh1_max)));
		
		//std::cout <<"k1, cos = " << k1 << ", "<<cosTh<<"\t" << Theta1_max << std::endl;
		
		const double factor_p1p2 = 1/(4*k0*Ecm);
		
		
		const vec3 e3(sinTh,0,cosTh);
		const vec3 e1(cosTh,0,-sinTh);
		const vec3 e2(0,1,0);	
		const vec3 vk1 = k1*e3;	
		const double deltaE = (k0-k1)*(k0+k1)/(Ek0+Ek1);
				
		return integrateAB([/*q_min,cosTh1_max,Theta1_max,NTh1,*/
							mp,k0,M2,vk1,Ek0,Ecm,Ek1,deltaE,St2,k1,Nphi,e1,e2,e3,factor_p1p2](double Theta1){
				
				const double cosTh1 = cos(Theta1);
				const double sinTh1 = sin(Theta1);
 
				const double q_part = Ecm - Ek1 + k1*cosTh1;
				const double q =  (Ecm*deltaE)/q_part;
				
				if(q_part <0)
					std::cout<< Ecm << ", " << Ek1 << ", k1 =" << k1 << std::endl;
				/*
				if(Theta1 < Theta1_max + M_PI/NTh1 && cosTh1_max < 1 && cosTh1_max > -1){
					std::cout<< "Theta1 = " << Theta1 << ", \nTheta1_max = " << Theta1_max << ", cos = " << cosTh1_max << 
						", \neps = " << q_min << ", \nq = " << q << std::endl;
				}
				*/
				const double phase_density = 
								( k1*k1*Ecm*deltaE)*sinTh1/(8*q_part*q_part*phase_2pi4*Ek1);
				
				return integrateAB([mp,M2,St2,Ek1,vk1,e1,e2,e3,sinTh1,cosTh1,q,factor_p1p2,phase_density](double phi1){
						
						const vec3 vq = q*(cosTh1*e3+sinTh1*sin(phi1)*e2+sinTh1*cos(phi1)*e1);

						PhaseState3 St3( vec4(Ek1,vk1),vec4(-vk1-vq,mp),vec4(q,vq));
						
							/*
							#define COS(a,b) (a.vecPart()*b.vecPart()/a.vecPart().norm()/b.vecPart().norm())
							if(phi1<2 && phi1 > 1 && Theta1>1 && Theta1 < 2 && cosTh>0.4 & cosTh <0.6 && k1 > 0.5*k0){
								std::cout << "q = " << q << ", vq = " << vq.norm() << std::endl;
								std::cout << "phi ="<<phi1<<",theta1 = "<< Theta1 <<", cosTh = "<<cosTh<<", k1 = "<< k1 <<std::endl;
								std::cout << cosTh << "vs" << COS(St3.K,St2.K) << std::endl;
								std::cout << cosTh1 << "vs" << COS(St3.K,St3.Q) << std::endl;
								std::cout << "St2 = "<<St2<<std::endl;
								std::cout << "St3 = "<<St3<<std::endl;
								std::cout << "Ecm = " << Ecm << std::endl;
								std::cout << "k1^2,p1^2,q^2 = " << 
												St3.K.quad() << ", " <<
												St3.P.quad() << ", " << 
												St3.Q.quad() << std::endl;
								std::cout <<(St3.K+St3.Q+St3.P) <<std::endl;
								
								std::cout <<"Ecm - K'-P'-Q"  << Ecm - (St3.K+St3.Q+St3.P).t << std::endl;
								std::cout << std::endl;
							}*/

						return M2(St3,St2)*phase_density*factor_p1p2;
						
					},0,2*M_PI,Nphi);
			},Theta1_max,M_PI,NTh1);
		};
	
}


extern inline dsigma_d3k1_type dsigma_d3k1_Q(double mk,double mp,double k0,double q_min,MatrixElementType23 M2,
								double Nphi ,double  NTh1){ 
	const double Ek0 = E(mk,k0);
	const double Ep0 = E(mp,k0);
	
	const double Ecm = Ek0+Ep0;
	
	
	const PhaseState2 St2(vec4(Ek0,0,0,k0),vec4(Ep0,0,0,-k0));
	
	const double k1p = ((1-q_min/Ecm)*sqrt(k0*k0+q_min*q_min+2*q_min*mk*mk/Ecm-2*q_min*Ek0) 
							- q_min*(Ek0-q_min)/Ecm)/(1-2*q_min/Ecm);
	
	const double k1m = ((1-q_min/Ecm)*sqrt(k0*k0+q_min*q_min+2*q_min*mk*mk/Ecm-2*q_min*Ek0) 
							+ q_min*(Ek0-q_min)/Ecm)/(1-2*q_min/Ecm);
	
	//std::cout<<"k1m = " << k1m << ",k1p = " << k1p << std::endl;
	//std::cout<<"cos_m = " << ( (Ecm-q_min)*(k0*k0-k1m*k1m)/(Ek0+E(mk,k1m)) - q_min*Ep0)/(q_min*k1m) << std::endl;
	//std::cout<<"cos_p = " << ( (Ecm-q_min)*(k0*k0-k1p*k1p)/(Ek0+E(mk,k1p)) - q_min*Ep0)/(q_min*k1p) << std::endl;

	return [q_min,mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,Nphi,NTh1,k1p,k1m](double k1,double cosTh){
		
		const double Ek1 = E(mk,k1);
		const double cosTh1_max = ( (Ecm-q_min)*(k0*k0-k1*k1)/(Ek0+Ek1) - q_min*Ep0)/(q_min*k1);
		
		const double sinTh = sqrt(1-cosTh*cosTh);
		
		const double Theta1_max = acos(MAX(-1,MIN(1,cosTh1_max)));
		
		//std::cout <<"k1, cos = " << k1 << ", "<<cosTh<<"\t" << Theta1_max << std::endl;
		
		const double factor_p1p2 = 1/(4*k0*Ecm);
		
		
		const vec3 e3(sinTh,0,cosTh);
		const vec3 e1(cosTh,0,-sinTh);
		const vec3 e2(0,1,0);	
		const vec3 vk1 = k1*e3;	
		const double deltaE = (k0-k1)*(k0+k1)/(Ek0+Ek1);
				
		return integrateAB5([/*q_min,cosTh1_max,Theta1_max,NTh1,*/
							mp,k0,M2,vk1,Ek0,Ecm,Ek1,deltaE,St2,k1,Nphi,e1,e2,e3,factor_p1p2](double Theta1){
				
				const double cosTh1 = cos(Theta1);
				const double sinTh1 = sin(Theta1);
 
				const double q_part = Ecm - Ek1 + k1*cosTh1;
				const double q =  (Ecm*deltaE)/q_part;
				
				if(q_part <0)
					std::cout<< Ecm << ", " << Ek1 << ", k1 =" << k1 << std::endl;
				/*
				if(Theta1 < Theta1_max + M_PI/NTh1 && cosTh1_max < 1 && cosTh1_max > -1){
					std::cout<< "Theta1 = " << Theta1 << ", \nTheta1_max = " << Theta1_max << ", cos = " << cosTh1_max << 
						", \neps = " << q_min << ", \nq = " << q << std::endl;
				}
				*/
				const double phase_density = 
								Ecm*deltaE*sinTh1/(8*q_part*q_part*phase_2pi5*Ek1);
				
				return integrateAB5([mp,M2,St2,Ek1,vk1,e1,e2,e3,sinTh1,cosTh1,q,factor_p1p2,phase_density](double phi1){
						
						const vec3 vq = q*(cosTh1*e3+sinTh1*sin(phi1)*e2+sinTh1*cos(phi1)*e1);

						PhaseState3 St3( vec4(Ek1,vk1),vec4(-vk1-vq,mp),vec4(q,vq));

						return M2(St3,St2)*phase_density*factor_p1p2;
						
					},0,2*M_PI,Nphi);
			},Theta1_max,M_PI,NTh1);
		};
	
}

extern inline double sigmaS(dsigma_dk1_dcosTh_type dsigma,double k,double k_esc,double k_tr,int Nk,int Nth){
	if(k_tr+k_esc > k){
		std::cout << "hard kinematics, use another function "<< std::endl;
		return 0;
	}
		
	
	
	
	if(k_esc > k_tr){
		double k_min = k_esc- k_tr;
		return integrateAB5([dsigma,Nth](double k1){
						return integrateAB5([dsigma,k1](double Th){
							return dsigma(k1,cos(Th))*sin(Th);
						},0,M_PI,Nth);
					},
			0,k_min,1 + (int)(Nk*k_min/2/k_esc)) + 

			integrateRightLogAB5([dsigma,Nth,k_esc,k_tr](double k1){
			double cosTh = (k_esc*k_esc-k1*k1-k_tr*k_tr)/(2*k1*k_tr);
			if(cosTh > 1 || cosTh < -1)
				std::cout << "cosTh = " << cosTh << " , k1 = " << k1 << std::endl;
			
				return integrateAB5([dsigma,k1](double Th){
					return dsigma(k1,cos(Th))*sin(Th);
					},acos(cosTh),M_PI,Nth);
				},k_min,k_esc+k_tr,k-k_esc-k_tr,Nk);
	}
	else{
		return integrateRightLogAB5([dsigma,Nth,k_esc,k_tr](double k1){
			double cosTh = (k_esc*k_esc-k1*k1-k_tr*k_tr)/(2*k1*k_tr);
			if(cosTh > 1 || cosTh < -1)
				std::cout << "cosTh = " << cosTh << " , k1 = " << k1 << std::endl;
			return integrateAB5([dsigma,k1](double Th){
				return dsigma(k1,cos(Th))*sin(Th);
				},acos(cosTh),M_PI,Nth);
			},k_tr-k_esc,k_esc+k_tr,k-k_esc-k_tr,Nk);
	}
	
	
}	
extern inline double sigmaC(dsigma_dk1_dcosTh_type dsigma,dsigma_dk1_dcosTh_type condition,double k,int Nk,int Nth){

	return integrateRightLogAB5([dsigma,Nth,condition](double k1){
		return integrateAB5([dsigma,k1,condition](double Th){
			return dsigma(k1,cos(Th))*sin(Th)*condition(k1,cos(Th));
			},0,M_PI,Nth);},  
		0,0.999*k,0.001*k,Nk);
}

extern inline double W(double x){
	const double y = x*x;
	if(y < 1./128){
		return 2./3.*y*(1. + 3./5*y*(
						1. + 5./7*y*(
						1. + 7./9*y*(
						1. + 9./11*y*(
						1. + 11./13*y*(
						1. + 13./15*y*(
						1. + 15./17*y)))))));
	}
	else
		return log((1.+x*(2.+x))/(1.-y))/x-2.;
}



extern inline double sigmaC(double mp,double mk,double v_ls,double v_esc,
							MatrixElementType22 M22,MatrixElementType23 M23,
							int Nphi = 6,int NTh1 = 6,int NTh = 10,int Nk = 40,
							double delta = 0.01){
							
	double k0 = v_ls*mp*mk/(mp+mk);
	double kt = v_ls*mk*mk/(mp+mk);
	double ke = v_esc*mk;
	
	double Ep = E(mp,k0);
	double Ek = E(mk,k0);
	
	double eps = delta*(Ep+Ek)*(k0*k0/(Ek+mk))/(Ek+Ep-mk);
	
	//std::cout << "k0 = " << k0 << std::endl;
	//std::cout << "q_min = " << eps << std::endl;
		
	double Ecm = Ep + Ek;
	
	const double k1p = ((1-eps/Ecm)*sqrt(k0*k0+eps*eps+2*eps*mk*mk/Ecm-2*eps*Ek) 
							- eps*(Ek-eps)/Ecm)/(1-2*eps/Ecm);
	
	const double k1m = ((1-eps/Ecm)*sqrt(k0*k0+eps*eps+2*eps*mk*mk/Ecm-2*eps*Ek) 
							+ eps*(Ek-eps)/Ecm)/(1-2*eps/Ecm);
							
	//auto dsigmaN = dsigma_d3k1_Q(mk,mp,k0,eps,M23,Nphi,NTh1);
	
	if(kt*(1+1e-8) >= ke + k0)
		return 0;
	if(kt+ke <= k0){
		//std::cout << "easy kinematics (kt+ke < k1p)" <<std::endl;
		auto dsigmaN = dsigma_d3k1(mk,mp,k0,M23,Nphi,NTh1);
		/*
		{
			auto dsigmaN_test = dsigma_dk1_dcosTh(mk,mp,k0,M23,Nphi,NTh1);
			std::cout << "easy kinematics (kt+ke < k1p), easy calculation" <<std::endl;
			std::cout << sigmaS(dsigmaN_test, k0,ke,kt,Nk,NTh)<< std::endl;
		}
		*/
		return integrateAB5([kt,dsigmaN,NTh](double ke_s){
				return integrateAB5([ke_s,kt,dsigmaN](double ThetaVe){
						double kx = ke_s*sin(ThetaVe);
						double kz = - kt + ke_s*cos(ThetaVe);
						
						double k1 = sqrt(kx*kx+kz*kz);
						double cosTh = 0;
						if(kz > 0)
							cosTh = 1./sqrt(1.+(kx*kx)/(kz*kz));
						else
							cosTh = -1./sqrt(1.+(kx*kx)/(kz*kz));
							
						return dsigmaN(k1,cosTh)*sin(ThetaVe)*2*M_PI*ke_s*ke_s;
					},0,M_PI,NTh);
			},0,ke,Nk);
	}
	else{
		return 0;
		//std::cout << "hard kinematics (kt+ke < k1p)" <<std::endl;
		double cosTh_cr = (ke*ke - kt*kt-k0*k0)/(2*k0*kt);
		//std::cout << "cosTh_cr = " << cosTh_cr <<std::endl;
		
		double Th_cr = acos(cosTh_cr);
		
		double alpha = 1./137.036;

		
		const vec4 P0 = vec4(vec3(0,0,-k0),mp);
		const vec4 K0 = vec4(vec3(0,0,k0),mk);
		
		//double reg_sigma = integrateAB5([Ecm,&P0,&K0,&M22,alpha,delta,mp,mk,NTh,k0](double Th){
		double reg_sigma = integrateAB5([Ecm,&P0,&K0,&M22,alpha,delta,mp,mk,NTh,k0](double Th){
				
				const vec4 P1 = vec4(vec3(-k0*sin(Th),0,-k0*cos(Th)),mp);
				const vec4 K1 = vec4(vec3(k0*sin(Th),0,k0*cos(Th)),mk);
				
				double M2 =  M22(PhaseState2(K1,P1),PhaseState2(K0,P0));
				
				double sc = P1*P0;
				double x = k0*sqrt((1.0-cos(Th))*(sc+mp*mp))/sc;
				double factor = 0.0+alpha/M_PI*log(delta)*W(x);
				
				return factor*M2/(32*M_PI*Ecm*Ecm)*sin(Th);
			}, Th_cr,M_PI,NTh);  
		
		//std::cout << "reg sigma = " <<reg_sigma <<std::endl;
		
		auto dsigmaN = dsigma_d3k1_Q(mk,mp,k0,eps,M23,Nphi,NTh1);
		
		//std::cout <<"eps_ = k0 - k1m = " <<k0-k1m << ", k1m-k1p = " << k1m-k1p <<std::endl; 
		
		double sigma_inter = integrateRightLogAB5([dsigmaN,ke,kt,NTh](double k1){
					double cosTh_min = (ke*ke-k1*k1-kt*kt)/(2*k1*kt);
					return integrateAB5([k1,dsigmaN](double Th){
							return sin(Th)*dsigmaN(k1,cos(Th))*k1*k1*2*M_PI;
						},acos(cosTh_min),M_PI,NTh);
			},k1p,k1m,k0-k1m,1 + (int)(Nk/8.0));
		
		//std::cout << "sigma intermediate from kp to km = " << sigma_inter << std::endl;
		
		double sigma_gamma = 0.0;
		double k_min = kt - ke;
		if(ke > kt){
			//std::cout << "ke > kt" <<std::endl;
			k_min = ke - kt;
			sigma_gamma += integrateAB5([dsigmaN,NTh](double k1){
					return k1*k1*2*M_PI*integrateAB5([dsigmaN,k1](double Th){
							return sin(Th)*dsigmaN(k1,cos(Th));
						},0,M_PI,NTh);
				},0,ke-kt,1 + (int)((ke-kt)/k0*Nk));
		}
		/*
		sigma_gamma +=  integrateAB5([k_min,k1p,NTh,ke,kt,dsigmaN](double tech_th_k1){
				double k1 = (k_min+k1p)/2 + (k1p-k_min)/2*cos(tech_th_k1);
				
				double cosTh_min = cosTh_min = (ke*ke-k1*k1-kt*kt)/(2*k1*kt);
				
				return 2*M_PI*k1*k1*(k1p-k_min)/2*sin(tech_th_k1)*
					integrateAB5([k1,dsigmaN](double Th){
						return sin(Th)*dsigmaN(k1,cos(Th));
					},acos(cosTh_min),M_PI,NTh);
			},0,M_PI,Nk);*/
			
		sigma_gamma +=  integrateRightLogAB5([NTh,ke,kt,dsigmaN](double k1){
				
				double cosTh_min = cosTh_min = (ke*ke-k1*k1-kt*kt)/(2*k1*kt);
					
				return 2*M_PI*k1*k1*
					integrateAB5([k1,dsigmaN](double Th){
						return sin(Th)*dsigmaN(k1,cos(Th));
					},acos(cosTh_min),M_PI,NTh);
			},k_min,k1p,k0-k1p,Nk);
			
		//std::cout << "sigma main = " << sigma_gamma << std::endl;
		return sigma_gamma + reg_sigma + sigma_inter;
	}
	
	
	
}

extern inline double sigmaC(double mp,double mk,double v_ls,double v_esc,MatrixElementType23 M23,
							int Nphi = 6,int NTh1 = 6,int NTh = 10,int Nk = 40){
							
	double k0 = v_ls*mp*mk/(mp+mk);
	double kt = v_ls*mk*mk/(mp+mk);
	double ke = v_esc*mk;
	
	double Ep = E(mp,k0);
	double Ek = E(mk,k0);

		
	double Ecm = Ep + Ek;
							
	
	if(kt*(1+1e-8) >= ke + k0)
		return 0;
	if(kt+ke <= k0){
		auto dsigmaN = dsigma_d3k1(mk,mp,k0,M23,Nphi,NTh1);
		return integrateAB5([kt,dsigmaN,NTh](double ke_s){
				return integrateAB5([ke_s,kt,dsigmaN](double ThetaVe){
						double kx = ke_s*sin(ThetaVe);
						double kz = - kt + ke_s*cos(ThetaVe);
						
						double k1 = sqrt(kx*kx+kz*kz);
						double cosTh = 0;
						if(kz > 0)
							cosTh = 1./sqrt(1.+(kx*kx)/(kz*kz));
						else
							cosTh = -1./sqrt(1.+(kx*kx)/(kz*kz));
							
						return dsigmaN(k1,cosTh)*sin(ThetaVe)*2*M_PI*ke_s*ke_s;
					},0,M_PI,NTh);
			},0,ke,Nk);
	}
	else{
		return 0;
	}
}

const double MX = RAND_MAX; 
const double FcE = 1e-10;
const double Fc = (1.0 - FcE);
extern inline double weight_xi(double xi,double xi_min){
	return sqrt(2/M_PI)*xi*exp(-xi_min*xi_min/2); 
}
extern inline double random_xi(double xi_min){
	const double dT = exp(-xi_min*xi_min/2);
	const double t = std::rand()/MX*Fc+FcE;
	return sqrt(-2*log(dT*t));
}
extern inline double random_cos(){
	return 2*std::rand()/MX-1;
}

extern inline double sigma_facotor_capture(double v,double u0,double cosTh_0, double cosTh1,int type= 0){
	double x = (1.0-cosTh1)/2;
	double y = (1.0+cosTh1)/2;
	double t = (v/u0)*(v/u0);
	if(type == 0)
		return x;
	else if(type == 1)
		return x*(1.0-y*cosTh_0)*t;
	else {
		return x*(1.0-y*
					(1.5*cosTh_0+0.75*cosTh_0*cosTh_0-0.25-
					y*(1.5*cosTh_0*cosTh_0-0.5)))*t*t;
	}
}

extern inline double sigma_facotor_esc(double v,double u0,double cosTh_0, double cosTh1,int type= 0){
	double x = (1.0-cosTh1)/2;
	double y = (1.0+cosTh1)/2;
	double t = (v/u0)*(v/u0);
	if(type == 0)
		return y;
	else if(type == 1)
		return y*(1.0+x*cosTh_0)*t;
	else {
		return y*(1.0+x*
					(1.5*cosTh_0-0.75*cosTh_0*cosTh_0+0.25+
					x*(1.5*cosTh_0*cosTh_0-0.5)))*t*t;
	}
}
 
enum EscapeOrCapture{
	CAPTURE,ESCAPE
};

double sigmaTfacor(double mp,double mk,double v,double vesc,double u0,double wT,
					EscapeOrCapture esc,int type,size_t N){
	
	double wTmin = (mp+mk)/(2*mp)*(std::abs((mp-mk)/(mp+mk))*v - vesc);
	
	
	
	if(wTmin <0)
		wTmin = 0;
	
	
	double rmin = wTmin/wT;
	if(wT <=0)
		rmin = 0;
	//std::cout << "wr = " << rmin <<std::endl;
	//std::cout << std::string("rmin = ") + std::to_string(rmin) + "\n";
	if(rmin > 10)
		return 0.0;
	
	double sum = 0.0;
	vec3 V(0,0,v);
	#pragma omp parallel reduction(+:sum)
	{
		double sumt = 0.0;
		size_t threadnum = omp_get_num_threads();
		size_t Nt = N/threadnum + 1;
		
		//std::cout << threadnum << std::endl;
		std::srand(int(time(NULL)) * (omp_get_thread_num()+1));
		//#pragma omp parallel for reduction(+:sum)
		for(size_t i=0;i<Nt;i++){
			double cosR = random_cos();
			double r = random_xi(rmin);
			double w = r*wT;
			vec3 W(0,w*sqrt(1-cosR*cosR),w*cosR);
			vec3 Vt = (mp*W+mk*V)/(mp+mk);
			vec3 Vc = (mp/(mp+mk))*(V-W);
			double vt = Vt.norm();
			double vc = Vc.norm();
			
			if(vc != 0){
				if(vesc+vt <= vc){
					if(esc == ESCAPE){
						sumt += weight_xi(r,rmin)*pow(vc/u0,type*2);
						//std::cout << "easy ESCAPE" <<std::endl;
					}
				}
				else if(vc+vt<=vesc){
					if(esc == CAPTURE)
						sumt += weight_xi(r,rmin)*pow(vc/u0,type*2);
						//std::cout << "easy CAPTURE" <<std::endl;
				}
				else if(vesc+vc <= vt){
					if(esc == ESCAPE)
						sumt += weight_xi(r,rmin)*pow(vc/u0,type*2);
				}
				else if(vt >= 0){
				
					double cosTh_0 = -(Vt*Vc)/(vt*vc);
					double cosTh_1 = (vt*vt+vc*vc-vesc*vesc)/(2*vc*vt);
					
					if(-1.0 <= cosTh_1 &&cosTh_1 <= 1.0){
						if(esc == CAPTURE){
							sumt += weight_xi(r,rmin)*sigma_facotor_capture(vc,u0,cosTh_0,cosTh_1,type);
						}
						else{
							sumt += weight_xi(r,rmin)*sigma_facotor_esc(vc,u0,cosTh_0,cosTh_1,type);
						}
					}
				}
			}
			//sumt +=  weight_xi(r,rmin)*w*cosR*cosR;
			/*
			if(1){
				std::cout << std::to_string(omp_get_thread_num()) + ", i= " + std::to_string(i) + 
					": " + std::to_string(sum) +  "\nr = " + std::to_string(r) + "\n\n";
			}*/
		}
		sum += sumt/(Nt*threadnum);
	}
	return sum;
}

#endif

