#ifndef MATRIX_ELEMENT_H
#define MATRIX_ELEMENT_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include "integrator.hpp" 
#include "phase4.hpp"
#include <functional>
#include "cross_section3.hpp"

const double e2 = (4*M_PI/137.036);

extern inline double q_coeff(vec4 p1,vec4 p2,vec4 q){
	vec4 dp = p2-p1;
	double p1q = p1*q;
	double p2q = p2*q;
	return (( (dp*q)*p1 - dp*p1q )/(p1q*p2q) ).quad();
}

MatrixElementType22 ScalarElement(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	return [a1,a2,b1,b2,mp,mk](PhaseState2 out, PhaseState2 in){
		return 4*((a2*a2+b2*b2)*out.K0*in.K0 + (a2*a2-b2*b2)*mk*mk)*
				( (a1*a1+b1*b1)*out.P0*in.P0 + (a1*a1-b1*b1)*mp*mp );
	};

}

MatrixElementType23 ScalarElementQ(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	return [a1,a2,b1,b2,mp,mk](PhaseState3 out, PhaseState2 in){
		double M2hi =  4*((a2*a2+b2*b2)*out.K1*in.K0 + (a2*a2-b2*b2)*mk*mk);
		
		double M2ee = ((a1*a1+b1*b1)*out.P1*in.P0 + (a1*a1-b1*b1)*mp*mp);
		
		double p1q = out.P1*out.Q;
		double pq = in.P0*out.Q;
		
		double dM2ee = (a1*a1+b1*b1)*(in.P0-out.P1)*out.Q;
		double qcoeff = -(out.P1/p1q - in.P0/pq).quad();
		
		std::cout<< qcoeff << " vs " << q_coeff(in.P0,out.P1,out.Q);
		
		double part3 =  (a1*a1+b1*b1)*(pq-p1q)*(1/p1q-1/pq);
		
		return M2hi*(qcoeff*(M2ee+dM2ee)+part3)*e2;
	};

}

MatrixElementType22 VectorElement(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	double mltp = (a1*a1+b1*b1)*(a2*a2+b2*b2);
	double mlt2 = a1*a1*a2*a2-b1*b1*b2*b2;
	double mlt4 = a1*a2*b1*b2;
	return [mltp,mlt2,mlt4,mp,mk](PhaseState2 out, PhaseState2 in){
		
		return 4*(mltp*2*((in.K0*out.P0)*(out.K0*in.P0) +(in.K0*in.P0)*(out.K0*out.P0) ) + 
				mltp*4*mk*mk*mp*mp - 4*mlt2*(mp*mp*(out.K0*in.K0) + mk*mk*(out.P0*in.P0)) + 
				4*mlt4*((in.K0*in.P0)*(out.K0*out.P0) - (in.K0*out.P0)*(out.K0*in.P0) ));
	};
}

MatrixElementType23 VectorElementQ(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	double mltp = (a1*a1+b1*b1)*(a2*a2+b2*b2);
	double mlt2 = a1*a1*a2*a2-b1*b1*b2*b2;
	double mlt4 = a1*a2*b1*b2;
	return [mltp,mlt2,mlt4,mp,mk](PhaseState3 out, PhaseState2 in){	
		double p1q = out.P1*out.Q;
		double pq = in.P0*out.Q;
		double qcoeff = -(out.P1/p1q - in.P0/pq).quad();
		return 4*qcoeff*(mltp*2*((in.K0*out.P1)*(out.K1*in.P0) +(in.K0*in.P0)*(out.K1*out.P1) ) + 
				mltp*4*mk*mk*mp*mp - 4*mlt2*(mp*mp*(out.K1*in.K0) + mk*mk*(out.P1*in.P0)) + 
				4*mlt4*((in.K0*in.P0)*(out.K1*out.P1) - (in.K0*out.P1)*(out.K1*in.P0) ));
	};

}

#endif
