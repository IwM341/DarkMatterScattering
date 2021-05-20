#ifndef CROSSSEC_H
#define CROSSSEC_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include "integrator.hpp" 
#include "phase4.hpp"
#include <functional>


#define  phase_2pi1   6.2831853071795864770
#define  phase_2pi2   39.478417604357434476
#define  phase_2pi3   248.05021344239856142
#define  phase_2pi4   1558.5454565440389959
#define  phase_2pi5   9792.6299131290065050

inline double fase_volume_density(double k1,double Ecm,double Ek0,double Ek1,double cosTh1){
	double q_part = Ecm - Ek1 - k1*cosTh1;
	return ( k1* k1*Ecm*(Ek0-Ek1))/(8*q_part*q_part*phase_2pi4*Ek1);
}

//pase volume parameters: dk' d cosTheta d cosTheta1 d phi1

//integration through: d cosTheta1 d phi1


class PhaseState3{
public:
	vec4 K1,P1,Q;
	PhaseState3(vec4 K1,vec4 P1,vec4 Q):K1(K1),P1(P1),Q(Q){}
	PhaseState3(){}
	

};
class PhaseState2{
public:
	vec4 K0,P0;
	PhaseState2(vec4 K0,vec4 P0):K0(K0),P0(P0){}
	PhaseState2(){}
	

};



std::ostream& operator<<(std::ostream& os, const PhaseState3& V3){
	os << "k1 = " << V3.K1.toString() << ", p1 = "<<V3.P1.toString()<<",q = " <<V3.Q.toString();
	return os;
}

std::ostream& operator<<(std::ostream& os, const PhaseState2& V2){
	os << "k0 = " << V2.K0.toString() << ", p0 = "<<V2.P0.toString();
	return os;
}

typedef std::function<double(PhaseState3,PhaseState2)> MatrixElementType;

typedef std::function<double(double k1,double cosTh)> dsigma_dk1_dcosTh_type;


extern inline dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh(double mk,double mp,double k0,MatrixElementType M2,
								double Nphi ,double  NTh1){ 
	double Ek0 = E(mk,k0);
	double Ep0 = E(mp,k0);
	
	double Ecm = Ek0+Ep0;
	
	
	PhaseState2 St2(vec4(Ek0,0,0,k0),vec4(Ep0,0,0,-k0));
	
	return [mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,Nphi,NTh1](double k1,double cosTh){
		return integrateAB([mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,k1,cosTh,Nphi,NTh1](double Theta1){
				return integrateAB([mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,k1,cosTh,Theta1,Nphi,NTh1](double phi1){
						
						
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
							
							
							St3.K1 = vec4(Ek1,vk1);
							St3.Q = vec4(q,vq);
							St3.P1 = vec4(-vk1-vq,mp);
							
							/*
							#define COS(a,b) (a.vecPart()*b.vecPart()/a.vecPart().norm()/b.vecPart().norm())
							if(phi1<2 && phi1 > 1 && Theta1>1 && Theta1 < 2 && cosTh>0.4 & cosTh <0.6 && k1 > 0.5*k0){
								std::cout << "q = " << q << ", vq = " << vq.norm() << std::endl;
								std::cout << "phi ="<<phi1<<",theta1 = "<< Theta1 <<", cosTh = "<<cosTh<<", k1 = "<< k1 <<std::endl;
								std::cout << cosTh << "vs" << COS(St3.K1,St2.K0) << std::endl;
								std::cout << cosTh1 << "vs" << COS(St3.K1,St3.Q) << std::endl;
								std::cout << "St2 = "<<St2<<std::endl;
								std::cout << "St3 = "<<St3<<std::endl;
								std::cout << "Ecm = " << Ecm << std::endl;
								std::cout << "k1^2,p1^2,q^2 = " << 
												St3.K1.quad() << ", " <<
												St3.P1.quad() << ", " << 
												St3.Q.quad() << std::endl;
								std::cout <<(St3.K1+St3.Q+St3.P1) <<std::endl;
								std::cout <<log ((St3.K1+St3.Q+St3.P1-St2.K0 - St2.P0).t) <<std::endl;
								std::cout <<"Ecm - K'-P'-Q"  << Ecm - (St3.K1+St3.Q+St3.P1).t << std::endl;
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


#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )

extern inline dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh_Q(double mk,double mp,double k0,double q_min,MatrixElementType M2,
								double Nphi ,double  NTh1){ 
	const double Ek0 = E(mk,k0);
	const double Ep0 = E(mp,k0);
	
	const double Ecm = Ek0+Ep0;
	
	
	const PhaseState2 St2(vec4(Ek0,0,0,k0),vec4(Ep0,0,0,-k0));
	
	const double k1p = ((1-q_min/Ecm)*sqrt(k0*k0+q_min*q_min+2*q_min*mk*mk/Ecm-2*q_min*Ek0) 
							- q_min*(Ek0-q_min)/Ecm)/(1-2*q_min/Ecm);
	
	const double k1m = ((1-q_min/Ecm)*sqrt(k0*k0+q_min*q_min+2*q_min*mk*mk/Ecm-2*q_min*Ek0) 
							+ q_min*(Ek0-q_min)/Ecm)/(1-2*q_min/Ecm);
	
	std::cout<<"k1m = " << k1m << ",k1p = " << k1p << std::endl;
	std::cout<<"cos_m = " << ( (Ecm-q_min)*(k0*k0-k1m*k1m)/(Ek0+E(mk,k1m)) - q_min*Ep0)/(q_min*k1m) << std::endl;
	std::cout<<"cos_p = " << ( (Ecm-q_min)*(k0*k0-k1p*k1p)/(Ek0+E(mk,k1p)) - q_min*Ep0)/(q_min*k1p) << std::endl;

	return [q_min,mk,mp,k0,M2,Ek0,Ep0,Ecm,St2,Nphi,NTh1,k1p,k1m](double k1,double cosTh){
		
		const double Ek1 = E(mk,k1);
		const double cosTh1_max = ( (Ecm-q_min)*(k0*k0-k1*k1)/(Ek0+Ek1) - q_min*Ep0)/(q_min*k1);
		
		const double sinTh = sqrt(1-cosTh*cosTh);
		
		const double Theta1_max = acos(MAX(-1,MIN(1,cosTh1_max)));
		
		std::cout <<"k1, cos = " << k1 << "\t" << cosTh1_max << std::endl;
		
		const double factor_p1p2 = 1/(4*k0*Ecm);
		
		
		const vec3 e3(sinTh,0,cosTh);
		const vec3 e1(cosTh,0,-sinTh);
		const vec3 e2(0,1,0);	
		const vec3 vk1 = k1*e3;	
		const double deltaE = (k0-k1)*(k0+k1)/(Ek0+Ek1);
		std::cout<<deltaE<<" vs " << Ek0-Ek1 <<std::endl;
				
		return integrateAB5([mp,k0,M2,vk1,Ek0,Ecm,Ek1,deltaE,St2,k1,Nphi,e1,e2,e3,factor_p1p2](double Theta1){
				
				const double cosTh1 = cos(Theta1);
				const double sinTh1 = sin(Theta1);
 
				const double q_part = Ecm - Ek1 + k1*cosTh1;
				const double q =  (Ecm*deltaE)/q_part;
				
				if(q_part <0)
							std::cout<< Ecm << ", " << Ek1 << ", k1 =" << k1 << std::endl;
				

				
				const double phase_density = 
								( k1*k1*Ecm*deltaE)*sinTh1/(8*q_part*q_part*phase_2pi4*Ek1);
				
				return integrateAB5([mp,M2,St2,Ek1,vk1,e1,e2,e3,sinTh1,cosTh1,q,factor_p1p2,phase_density](double phi1){
						
						
						const vec3 vq = q*(cosTh1*e3+sinTh1*sin(phi1)*e2+sinTh1*cos(phi1)*e1);
						
						
						

						PhaseState3 St3( vec4(Ek1,vk1),vec4(-vk1-vq,mp),vec4(q,vq));
						
							/*
							#define COS(a,b) (a.vecPart()*b.vecPart()/a.vecPart().norm()/b.vecPart().norm())
							if(phi1<2 && phi1 > 1 && Theta1>1 && Theta1 < 2 && cosTh>0.4 & cosTh <0.6 && k1 > 0.5*k0){
								std::cout << "q = " << q << ", vq = " << vq.norm() << std::endl;
								std::cout << "phi ="<<phi1<<",theta1 = "<< Theta1 <<", cosTh = "<<cosTh<<", k1 = "<< k1 <<std::endl;
								std::cout << cosTh << "vs" << COS(St3.K1,St2.K0) << std::endl;
								std::cout << cosTh1 << "vs" << COS(St3.K1,St3.Q) << std::endl;
								std::cout << "St2 = "<<St2<<std::endl;
								std::cout << "St3 = "<<St3<<std::endl;
								std::cout << "Ecm = " << Ecm << std::endl;
								std::cout << "k1^2,p1^2,q^2 = " << 
												St3.K1.quad() << ", " <<
												St3.P1.quad() << ", " << 
												St3.Q.quad() << std::endl;
								std::cout <<(St3.K1+St3.Q+St3.P1) <<std::endl;
								
								std::cout <<"Ecm - K'-P'-Q"  << Ecm - (St3.K1+St3.Q+St3.P1).t << std::endl;
								std::cout << std::endl;
							}*/

						return M2(St3,St2)*phase_density*factor_p1p2;
						
					},0,2*M_PI,Nphi);
			},0,M_PI,NTh1);
		};
	
}
#endif

