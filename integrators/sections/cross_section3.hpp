#ifndef CROSSSEC_H
#define CROSSSEC_H

#include "phase4.hpp"
#include "lambda_int.hpp"

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

typedef std::function<double(PhaseState3,PhaseState2)> MatrixElementType;

typedef std::function<double(double k1,double cosTh)> dsigma_dk1_dcosTh_type;

dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh(double mk,double mp,double k0,MatrixElementType M2,
								double Nphi = 100,double  NTh1 = 100);
								
dsigma_dk1_dcosTh_type dsigma_dk1_dcosTh_fast(double mk,double mp,double k0,MatrixElementType M2,
								double Nphi = 100,double  NTh1 = 100);


#endif
