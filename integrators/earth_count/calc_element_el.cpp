#include <iostream>
#include "csv_function.hpp"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include "utils.hpp"

#define SVAR(x) (std::string(#x) + std::string(" = ") + std::to_string(x))
#define PVAR(x) std::cout << std::string(#x) + std::string(" = ") + std::to_string(x) <<std::endl;

namespace std{
	string to_string(const string & S){
		return S;
	}
}

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


Function1<double> MassElement(std::string el,bool isElastic,int sigma_type,
					double Mmin,double Mmax,size_t N ,double q,double sigma_U,int isTemp,
					const std::string&path = ""){
	auto &EM = BodyModel::Instance();
	double Mres = ME.at(el);
	double MmaxO = Mmax;
	if(!isElastic){
		Mmax = Mres;
		Mres = (Mmax+Mmin)/2;
		
	}
	
	int Nmk = 1;
	if(isTemp){
		Nmk = 100000;
	}
	
	std::vector<double> mgrid = Mgrid(Mmin,Mmax,Mres,N,q,true);
	
	if(mgrid[mgrid.size()-1] < MmaxO-0.001){
		mgrid.push_back(MmaxO);
	}
	
	std::vector<double> CoeffGrid(mgrid.size());
	std::vector<double> VescGrid = grid(20,EM.VeMin(),EM.VeMax(),1,0.0);
	std::vector<double> Ugrid;
	if(isElastic)
		Ugrid = grid(300,EM.VeMin()/100,6*U0,4,0.001);
	else
		Ugrid = grid(100,U0/10,U0*6,1,0.0);
	
	//std::cout << mgrid <<std::endl;
	Function2<double> Sgm;
	for(size_t i=0;i<mgrid.size();i++){
		//std::cout << mgrid[i] <<std::endl;
		if(isElastic)
			Sgm = SigmaElastic(mgrid[i],el,Ugrid,VescGrid,isTemp,sigma_type,Nmk);
		else
			Sgm = SigmaInelastic(mgrid[i],el,Ugrid,VescGrid,sigma_type);
		Function1<double> Su = IntegrateSigma(Sgm,el);
		
		CoeffGrid[i] = C_ND(Su,sigma_U*U0);
	}
	
	
		
	return Function1<double>(mgrid,CoeffGrid);
	/*
	std::ofstream outf(path + el + "(" + Elasticness[isElastic] + ", " + 
								TempConsider[isTemp] + std::to_string(sigma_type) + ").dat");
	outf << "Mw\tCnd" << std::endl;
	outf << Function1<double>(mgrid,CoeffGrid) << std::endl;
	*/
}



std::string Filename(std::string el,bool isElastic,int sigma_type,
					double Mmin,double Mmax,size_t N ,double q,double sigma_U,int isTemp,
					const std::string&path = ""){
	
	std::string Elasticness[2]  = {"in","el"};
	std::string TempConsider[2]  = {"","T, "};				
	return path + el + "(" + Elasticness[isElastic] + ", " + 
								TempConsider[isTemp] + std::to_string(sigma_type) + ").dat";
}

Function1<double> Divide(Function1<double> IN,Function1<double> EL){
	Function1<double> DIV = IN;
	DIV.insert(EL);
	DIV = Function1<double>(DIV.getXgrid(),[&IN,&EL](double mass){return IN(mass)/EL(mass);});
	return DIV;
}

int main(void){
	const std::vector<std::string> Els({"Fe","Ni","Mg","Si","O","Ca","Al","S","Cr","Na","P","Mn","C","H","Li"});
	const double sigmaU = 1.174;
	for(int sigmaType =0;sigmaType<3;sigmaType++){
		std::vector<double> Mgrid = grid(1000,2,100,1,0);
		std::vector<double> P(Mgrid.size(),0.0);
		for(auto el : Els){
			
	
			std::cout << SVAR(el) << "\t" <<SVAR(sigmaType) << std::endl;
			
			auto EL = MassElement(el,1,sigmaType,2,100,100 ,1.6,sigmaU,0,std::string("resultE\\elastic\\") + std::to_string(sigmaType) +"\\");
			
			std::ofstream outEL(Filename(el,1,sigmaType,2,100,100 ,1.6,sigmaU,0, 
										std::string("resultE\\") + std::to_string(sigmaType) +"\\" ));
			outEL << "Mw\tCnd" << std::endl;
			outEL << EL << std::endl;
			
			double N = ME.at(el);
			P  = P + apply_function<double,double>(Mgrid,[N,&EL](double mk){
				return EL(mk)*N*(1+mk)*(1+mk)/(N+mk)/(N+mk)/mk;
			});
			
			/*
			auto DV = Divide(IN,EL);
			
			std::ofstream outDV(std::string("result\\divison\\") + std::to_string(sigmaType) +"\\" + 
									el + "(" + std::to_string(sigmaType) + ").dat");
			outDV << "Mw\tDiv" << std::endl;
			outDV << DV << std::endl;
			*/
			
		}
		
		std::ofstream out(std::string("resultE\\") + std::to_string(sigmaType) +"\\Combined.dat");
		out << "Mw\tP"<<std::endl;
		out << Function1<double>(Mgrid,P);
	}
}
/*
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
	
	auto F2 = MassElement(el,isElastic,sigma_type, Mmin, Mmax, N , q, sigma_U, isTemp );
	
	std::ofstream outf(Filename((el,isElastic,sigma_type, Mmin, Mmax, N , q, sigma_U, isTemp ));
	outf << "Mw\tCnd" << std::endl;
	outf << F2 << std::endl;
	
	
	
	
	return 0;
}
*/