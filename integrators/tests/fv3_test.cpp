#include "mc.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <phase4.hpp>
#include "integrators\lambda_int.hpp"

double deltaexp(double x,double sgm){
	return std::exp(-x*x/(2*sgm*sgm))/std::sqrt(2*M_PI*sgm*sgm);
}

double delta(double x,double sgm){
	if( x >= sgm/2 || x <= -sgm/2)
		return 0;
	else
		return 1.0/sgm;
}

class Volume3{
public:
	vec4 K,P,Q;
};

double rand_ab(double a,double b){
	return a + (std::rand()*(b-a))/RAND_MAX;
}

class v3_gen:public AreaGenerator<Volume3>{
private:
	double Ecm;
	double mp;
	double mk;
	
public:
	v3_gen(double Ecm,double mp,double mk):
	Ecm(Ecm),mp(mp),mk(mk){ srand(time(0));}
	
	virtual Volume3 RandomCoords() const{
		vec3 k1(rand_ab(-Ecm,Ecm),rand_ab(-Ecm,Ecm),rand_ab(-Ecm,Ecm));
		
		double q = rand_ab(0,Ecm);
		
		vec3 q1 = vec3::PolarCos(q, rand_ab(-1,1),rand_ab(0,2*M_PI));
		
		Volume3 V;
		V.K = vec4(k1,mk);
		V.Q = vec4(q,q1);
		V.P = vec4(-(q1 + k1),mp);
		return V;
	}
	virtual double AreaVolume() const{return 8*Ecm*Ecm*Ecm*Ecm*4*M_PI;}
	double density(Volume3 V) const{
		return pow(1.0/(2*M_PI),5)*V.Q.t/(8*V.K.t*V.P.t)*delta(V.K.t+V.P.t+V.Q.t-Ecm,Ecm/500);
	}
};

typedef  std::function<double(Volume3)> ftype;

inline double kdensity(double k,double Ecm,double Ek0,double mk){
	double dE = Ek0-E(mk,k);
	double Ep = Ecm-E(mk,k);
	
	return  (k*k*Ecm*dE)/(2*pow(2*M_PI,3) * E(mk,k)*(  Ep*Ep-k*k  )); 	
}

int main(void){
	
	vec3 k0(0,0,1);
	double mk = 2;
	double mp = 1;
	
	double Ecm = E(mk,k0.norm())+E(mp,k0.norm());
	
	double Ek0 = E(mk,k0.norm());
	
	std::cout << "true = " << 
		integrateAB([Ek0,mk,Ecm](double k){
			return kdensity(k,Ecm,Ek0,mk);
		},0,k0.norm(),200) << std::endl;

	
	MonteCarloIntegrator I;
	v3_gen gen(Ecm,mp,mk);
	std::cout<< I((ftype)[gen](Volume3 V)->double{return gen.density(V);},&gen,100000000) << std::endl;
	
	return 0;
}
