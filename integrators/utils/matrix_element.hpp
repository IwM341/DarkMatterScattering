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
	return -(( (dp*q)*p1 - dp*p1q )/(p1q*p2q) ).quad();
}



MatrixElementType22 MET0(double Ecm,double sigma0 = 1){
	double mult = sigma0*Ecm*Ecm*16*M_PI;
	return [mult](const PhaseState2 &out,const PhaseState2& in){
		return mult;
	};
}	
MatrixElementType23 MET0Q(double Ecm,double sigma0 = 1){
	double mult = sigma0*Ecm*Ecm*16*M_PI;
	return [mult](const PhaseState3 &out,const PhaseState2& in){
		return mult*q_coeff(in.P,out.P,out.Q)*e2;
	};
}

MatrixElementType22 MET1(double Ecm,double p0,int type =0,double sigma0 = 1 ){
	double mult = sigma0*Ecm*Ecm*4*M_PI/(p0*p0);
	if(type == 0){
		return [mult](const PhaseState2 &out,const PhaseState2 &in){
			//std::cout << "State["<< out << ", " << in << "]\n" <<std::endl; 
			return -mult*(out.P-in.P).quad();
		};
	}
	else{
		return [mult](const PhaseState2 &out,const PhaseState2 &in){
			return -mult*(out.K-in.K).quad();
		};
	}
	
}

MatrixElementType23 MET1Q(double Ecm,double p0,int type =0 ,double sigma0 = 1){
	double mult = sigma0*Ecm*Ecm*4*M_PI/(p0*p0);
	if(type == 0){
		return [mult](const PhaseState3 &out, PhaseState2 in){
			return -mult*(out.P-in.P).quad()*q_coeff(in.P,out.P,out.Q)*e2;
		};
	}
	else{
		return [mult](const PhaseState3 &out, PhaseState2 in){
			return -mult*(out.K-in.K).quad()*q_coeff(in.P,out.P,out.Q)*e2;
		};
	}
}
MatrixElementType22 MET2(double Ecm,double p0,double k0,double sigma0 = 1){
	double mult = sigma0*Ecm*Ecm*3*M_PI/(2*p0*p0*k0*k0);
	return [mult](const PhaseState2 &out, PhaseState2 in){
		return mult*(out.P-in.P).quad()*(out.K-in.K).quad();
	};
	
}
MatrixElementType23 MET2Q(double Ecm,double p0,double k0,double sigma0 = 1){
	double mult = sigma0*Ecm*Ecm*3*M_PI/(2*p0*p0*k0*k0);
	return [mult](const PhaseState3 &out, PhaseState2 in){
		return mult*(out.P-in.P).quad()*(out.K-in.K).quad()*
			q_coeff(in.P,out.P,out.Q)*e2;
	};
}

MatrixElementType22 ScalarElement(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	return [a1,a2,b1,b2,mp,mk](const PhaseState2 &out,const PhaseState2& in){
		return 4*((a2*a2+b2*b2)*out.K*in.K + (a2*a2-b2*b2)*mk*mk)*
				( (a1*a1+b1*b1)*out.P*in.P + (a1*a1-b1*b1)*mp*mp );
	};

}

MatrixElementType23 ScalarElementQ(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	return [a1,a2,b1,b2,mp,mk](const PhaseState3 &out,const PhaseState2& in){
		double M2hi =  4*((a2*a2+b2*b2)*out.K*in.K + (a2*a2-b2*b2)*mk*mk);
		
		double M2ee = ((a1*a1+b1*b1)*out.P*in.P + (a1*a1-b1*b1)*mp*mp);
		
		double p1q = out.P*out.Q;
		double pq = in.P*out.Q;
		
		double dM2ee = (a1*a1+b1*b1)*(in.P-out.P)*out.Q;
		double qcoeff = -(out.P/p1q - in.P/pq).quad();
		
		//std::cout<< qcoeff << " vs " << q_coeff(in.P,out.P,out.Q);
		
		double part3 =  (a1*a1+b1*b1)*(pq-p1q)*(1/p1q-1/pq);
		
		return M2hi*(qcoeff*(M2ee+dM2ee)+part3)*e2;
	};

}

MatrixElementType22 VectorElement(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	double mltp = (a1*a1+b1*b1)*(a2*a2+b2*b2);
	double mlt2 = a1*a1*a2*a2-b1*b1*b2*b2;
	double mlt4 = a1*a2*b1*b2;
	return [mltp,mlt2,mlt4,mp,mk](const PhaseState2 &out,const PhaseState2& in){
		
		return 4*(mltp*2*((in.K*out.P)*(out.K*in.P) +(in.K*in.P)*(out.K*out.P) ) + 
				mltp*4*mk*mk*mp*mp - 4*mlt2*(mp*mp*(out.K*in.K) + mk*mk*(out.P*in.P)) + 
				4*mlt4*((in.K*in.P)*(out.K*out.P) - (in.K*out.P)*(out.K*in.P) ));
	};
}

MatrixElementType23 VectorElementQ(double mp,double mk,double a1 = 1,double a2 = 1,double b1 = 0,double b2 = 0){
	
	double mltp = (a1*a1+b1*b1)*(a2*a2+b2*b2);
	double mlt2 = a1*a1*a2*a2-b1*b1*b2*b2;
	double mlt4 = a1*a2*b1*b2;
	return [mltp,mlt2,mlt4,mp,mk](const PhaseState3 &out,const PhaseState2& in){	
		double p1q = out.P*out.Q;
		double pq = in.P*out.Q;
		double qcoeff = -(out.P/p1q - in.P/pq).quad();
		return 4*qcoeff*(mltp*2*((in.K*out.P)*(out.K*in.P) +(in.K*in.P)*(out.K*out.P) ) + 
				mltp*4*mk*mk*mp*mp - 4*mlt2*(mp*mp*(out.K*in.K) + mk*mk*(out.P*in.P)) + 
				4*mlt4*((in.K*in.P)*(out.K*out.P) - (in.K*out.P)*(out.K*in.P) ));
	};

}

#endif
