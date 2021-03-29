#include "mc.hpp"
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <phase4.hpp>

double delta(double x,double sgm){
	return std::exp(-x*x/(2*sgm*sgm))/std::sqrt(2*M_PI*sgm*sgm);
}

double f(double x){
	return delta(x-0.5,0.01);
}

class Volume2{
public:
	vec4 K,P;
};

double rand_ab(double a,double b){
	return a + (std::rand()*(b-a))/RAND_MAX;
}

class v3_gen:public AreaGenerator<Volume2>{
private:
	double Ecm;
	double mp;
	double mk;
	
public:
	v3_gen(double Ecm,double mp,double mk):
	Ecm(Ecm),mp(mp),mk(mk){ srand(time(0));}
	
	virtual Volume2 RandomCoords() const{
		vec3 k1(rand_ab(-Ecm,Ecm),rand_ab(-Ecm,Ecm),rand_ab(-Ecm,Ecm));
		
		Volume2 V;
		V.K = vec4(k1,mk);
		V.P = vec4(-k1,mp);
		return V;
	}
	virtual double AreaVolume() const{return 8*Ecm*Ecm*Ecm;}
	double density(Volume2 V) const{
		return pow(1.0/(2*M_PI),2)/(V.K.t*V.P.t)*delta(V.K.t+V.P.t-Ecm,Ecm/100);
	}
};

typedef  std::function<double(Volume2)> ftype;

int main(void){
	
	vec3 k0(0,0,1);
	double mk = 1;
	double mp = 2;
	
	double Ecm = E(mk,k0.norm())+E(mp,k0.norm());
	
	std::cout << "true = " << 1/M_PI*k0.norm()/Ecm << std::endl;

	
	MonteCarloIntegrator I;
	v3_gen gen(Ecm,mp,mk);
	std::cout<< I((ftype)[gen](Volume2 V)->double{return gen.density(V);},&gen,10000000) << std::endl;
	
	return 0;
}
