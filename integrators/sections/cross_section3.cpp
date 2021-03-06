#include <cstdio>
#include <iostream>
#include <cmath>
#include "cross_section3.hpp"
#include "lambda_int.hpp" 

const double phase_2pi1 = 6.2831853071795864770; 
const double phase_2pi2 = 39.478417604357434476; 
const double phase_2pi3 = 248.05021344239856142; 
const double phase_2pi4 = 1558.5454565440389959; 
const double phase_2pi5 = 9792.6299131290065050; 

inline double fase_volume_density(double k1,double Ecm,double Ek0,double Ek1,double cosTh1){
	double q_part = Ecm - Ek1 - k1*cosTh1;
	return ( k1* k1*Ecm*(Ek0-Ek1))/(8*q_part*q_part*phase_2pi4*Ek1);
}

//pase volume parameters: dk' d cosTheta d cosTheta1 d phi1

//integration through: d cosTheta1 d phi1






std::ostream& operator<<(std::ostream& os, const PhaseState3& V3){
	os << "k1 = " << V3.K1.toString() << ", p1 = "<<V3.P1.toString()<<",q = " <<V3.Q.toString();
	return os;
}

std::ostream& operator<<(std::ostream& os, const PhaseState2& V2){
	os << "k0 = " << V2.K0.toString() << ", p0 = "<<V2.P0.toString();
	return os;
}



dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh(double mk,double mp,double k0,MatrixElementType M2,
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
						

						PhaseState3 St3;
						{
							double sinTh = sqrt(1-cosTh*cosTh);
							
							vec3 e3(sinTh,0,cosTh);
							vec3 e1(cosTh,0,-sinTh);
							vec3 e2(0,1,0);
							
							vec3 vk1 = k1*e3;
							
							
							double q =  (Ecm*(Ek0-Ek1))/(Ecm- Ek1 + k1*cosTh1);
							
							vec3 vq = q*(cosTh1*e3+sinTh1*sin(phi1)*e2+sinTh1*cos(phi1)*e1);
							
							
							St3.K1 = vec4(Ek1,vk1);
							St3.Q = vec4(q,vq);
							St3.P1 = vec4(-vk1-vq,mp);
							/*
							if(phi1<2 && phi1 > 1 && Theta1>1 && Theta1 < 2 && cosTh>0.4 & cosTh <0.6 && k1 > 0.5*k0){
								std::cout << "phi ="<<phi1<<",theta1 = "<< Theta1 <<", cosTh = "<<cosTh<<", k1 = "<< k1 <<std::endl;
								std::cout << "St2 = "<<St2<<std::endl;
								std::cout << "St3 = "<<St3<<std::endl;
							}*/
						}
						
						double q_part = Ecm - Ek1 + k1*cosTh1;
						double phase_density = 
								( k1*k1*Ecm*(Ek0-Ek1))*sinTh1/(8*q_part*q_part*phase_2pi4*Ek1);
						double factor_p1p2 = 1/(4*k0*Ecm);
						return M2(St3,St2)*phase_density*factor_p1p2;
						
					},0,2*M_PI,Nphi);
			},0,M_PI,NTh1);
		};
	
}



