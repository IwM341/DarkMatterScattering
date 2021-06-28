#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils.hpp"

#define SVAR(x) (std::string(#x) + std::string(" = ") + std::to_string(x))
#define PVAR(x) std::cout << std::string(#x) + std::string(" = ") + std::to_string(x) <<std::endl;


std::vector<double> Mgrid(double Mmin,double Mmax,double Mres,size_t N,double q = 1,bool eq = true){
	std::vector<double> M(N);
	
	if(Mres<Mmin){
		for(size_t i=0;i<N;i++)
			M[i] = Mmin+i*(Mmax-Mmin)/(N-1);
		return M;
	}
	
	double hmin = (Mres-Mmin)*
					pow( (1+pow((Mmax-Mres)/
								(Mres-Mmin),
							1/q))/(N-1),q);
	size_t N1 = 0.5+pow((Mres-Mmin)/hmin,1/q);
	
	M[N1] = Mres;
	
	if(eq){
		for(int i=1;i<N1+1;i++){
			M[N1-i] = Mres-hmin*pow(i,q);
		}
		for(int i=1;i<N-N1;i++){
			M[N1+i] = Mres+hmin*pow(i,q);
		}
		M[0] = Mmin;
		M[N-1] = Mmax;
	}
	else{
		for(int i=1;i<N1+1;i++){
			M[N1-i] = Mres-(Mres-Mmin)*pow((double)i/N1,q);
		}
		for(int i=1;i<N-N1;i++){
			M[N1+i] = Mres+(Mmax-Mres)*pow((double)i/(N-N1-1),q);
		}
	}
	return M;
}

std::string to_lower(const std::string &s) {
	std::string s1 = s;
    std::transform(s1.begin(), s1.end(), s1.begin(), 
                   [](unsigned char c){ return std::tolower(c); } // correct
                  );
    return s1;
}

template <class T1,class T2>
bool In(const T1 &x,const T2 &X){
	for(auto it = X.begin();it < X.end();it++){
		if(x == *it){
			return true;
		}
	}
	return false;
}

int main(int argc,char**argv){
	if(argc < 2){
		std::cout << "expect: [element name] [(isElastic)] [(sigma_type)] [(Mmin)] [(Mmax)] [(N)] [(q)] [(sigma_U)] [(Temp)]" << std::endl;
		return 0;
	}
	double Mmin = 2;
	double Mmax = 100;
	int sigma_type = 0;
	size_t N = 100;
	double q = 1;
	bool isElastic;
	double sigma_U = 1;
	int isTemp = 0;
	std::string el;
	
	el = argv[1];
	auto &EM = BodyModel::Instance();
	if(!EM.isExist(el) || ME.find(el) == ME.end()){
		std::cout << "not found element " << el << std::endl;
		return 0;
	}
	
	std::vector<std::string> args(argc - 1);
	for(size_t i = 1;i<argc;i++){
		args[i-1] = argv[i];
	}
	if(args.size() >= 2){
		if(In(to_lower(args[1]),std::vector<std::string>({"elastic","el","yes"})))
			isElastic = true;
		else if(In(to_lower(args[1]),std::vector<std::string>({"inelastic","in","no"})))
			isElastic = false;
		else{
			std::cout << "not understand " << "\"" << args[1] << "\"" << std::endl;
			return 0;
		}
		
		if(args.size() >= 3){
			sigma_type = std::stoi(args[2]);
			if(args.size() >= 5){
				Mmin = std::stod(args[3]);
				Mmax = std::stod(args[4]);
				if(args.size()>=6){
					N = std::stoi(args[5]);
					if(args.size()>=7){
						q = std::stod(args[6]);
						if(args.size()>=8){
							sigma_U = std::stod(args[7]);
							if(args.size()>=9){
								if(In(to_lower(args[8]),
										std::vector<std::string>({"t","yes","temp"})
									)){isTemp = 1;}
							}
						}
					}
				}
			}
		}
	}
	
	
	double Mres = ME.at(el);
	
	PVAR(Mmin);
	PVAR(Mres);
	PVAR(Mmax);
	PVAR(sigma_type);
	PVAR(N);
	PVAR(q);
	PVAR(isElastic);
	PVAR(isTemp);
	PVAR(sigma_U);
	
	if(!isElastic){
		Mmax = Mres;
		Mres = (Mmax+Mmin)/2;
		
	}
	
	std::vector<double> mgrid = Mgrid(Mmin,Mmax,Mres,N,q,true);
	std::vector<double> CoeffGrid(mgrid.size());
	std::vector<double> VescGrid = grid(20,EM.VeMin(),EM.VeMax(),1,0.0);
	std::vector<double> Ugrid;
	if(isElastic)
		Ugrid = grid(64,0.0,2*U0,6,0.01);
	else
		Ugrid = grid(64,U0/10,U0*3,1,0.0);
	
	//std::cout << mgrid <<std::endl;
	Function2<double> Sgm;
	for(size_t i=0;i<mgrid.size();i++){
		std::cout << mgrid[i] <<std::endl;
		if(isElastic)
			Sgm = SigmaElastic(mgrid[i],el,Ugrid,VescGrid,isTemp,sigma_type,10000);
		else
			Sgm = SigmaInelastic(mgrid[i],el,Ugrid,VescGrid,sigma_type);
		Function1<double> Su = IntegrateSigma(Sgm,el);
		CoeffGrid[i] = C_ND(Su,sigma_U*U0);
	}
	
	std::string Elasticness[2]  = {"in","el"};
	std::string TempConsider[2]  = {"","T, "};
	
	std::ofstream outf(el + "(" + Elasticness[isElastic] + ", " + 
								TempConsider[isTemp] + std::to_string(sigma_type) + ").dat");
	outf << "Mw\tCnd" << std::endl;
	outf << Function1<double>(mgrid,CoeffGrid) << std::endl;
	
	
	return 0;
}
