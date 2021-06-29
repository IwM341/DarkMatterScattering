#ifndef UTILS_H
#define UTILS_H
#include "csv_function.hpp"
#include "integrator.hpp"
#include "cross_section3.hpp"
#include "matrix_element.hpp"

#define U0 0.7667e-3



const std::map<std::string, double> ME({
			{"H1",1.0},{"He4",4.0},{"He3",3.0},{"C12",12.0},{"C13",13.0},
			{"N14",14.0},{"N15",15.0},{"O",16.0},{"O16",16.0},{"O17",17.0},{"O18",18.0},
			{"Ne",20.0},{"Na",23.0},{"Mg",24.0},{"Al",27.0},{"Si",28.0},{"P",31.0},
			{"S",32.0},{"Cl35",35.0},{"Cl36",36.0},{"Ar",40.0},{"K",39.0},{"Ca",40.0},
			{"Sc",45.0},{"Ti",48.0},{"V",51.0},{"Cr",52.0},
			{"Mn",55.0},{"Fe",56.0},{"Ni",59.0},{"Co",59.0},
			});

const std::map<std::string, int> QE({
			{"H1",1},{"He4",2},{"He3",2},{"C12",6},{"C13",6},
			{"N14",7},{"N15",7},{"O",8},{"O16",8},{"O17",8},{"O18",8},
			{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},
			{"S",16},{"Cl35",18},{"Cl36",18},{"Ar",18},{"K",19},{"Ca",20},
			{"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Ca",40.0},
			{"Mn",25},{"Fe",26},{"Ni",28}
			});

const std::vector<std::string> EarthElements({"Fe","O","Si","Mg","Ni","Ca","Al"}); 
	
class BodyModel{
	std::map<std::string, std::vector<double>> BM;
	double Vesc;
	double VescMax;
	BodyModel():BM(Function::CSVTable("earth_model.dat")),
		Vesc(3.7336e-5){
		
		if(BM.find("Vesc") == BM.end()){
			double vesc = Vesc;
			BM["Vesc"] = apply_function<double,double>
			(BM.at("phi"),[vesc](double x)->double{return sqrt(x)*vesc;});
		}
		VescMax = BM["Vesc"][0];
	}
	BodyModel(const BodyModel& root);
	BodyModel& operator = (const BodyModel&);
	
	
public:
	static BodyModel &Instance(){
		static BodyModel Inst;
		return Inst;
	}
	
	const std::vector<double> &operator[](const std::string & column) const{
		return BM.at(column);
	}
	bool isExist(const std::string & column) const{
		if(BM.find(column) == BM.end())
			return false;
		else 
			return true;
	}
	
	double VeMin() const{
		return Vesc;
	}
	double VeMax() const{
		return VescMax;
	}
	
};

Function1<double> IntegrateSigma(Function2<double> SigmaVVesc,const std::string &element){
	Function1<double> Sigma(SigmaVVesc.getXgrid());
	size_t N =  Sigma.size();
	
	const auto &PM = BodyModel::Instance();
	
	std::vector<double> US = SigmaVVesc.getXgrid();
	
	auto Vesc = PM["Vesc"];
	auto Xi = PM["Radius"];
	auto Rho = PM["RhoND"];
	auto Element = PM[element];
	
	
	double mass = ME.at(element);
	
	#pragma omp parallel for 
	for(size_t i=0;i<N;i++){
		double U = US[i];
		auto Vi = SigmaVVesc.atX(i);
		Sigma.getY(i) = 3.0/(U0*U0)*
		integrateAB2(Xi,
			apply_function<double,double>(Vesc,[U,&Vi](double vesc){
				double v = sqrt(vesc*vesc+U*U);
				return v*v*Vi(vesc);
				})*Xi*Xi*Element*Rho
		);
	}
	return Sigma;
}

Function2<double> SigmaInelastic(double mass,const std::string &element,
	const std::vector<double> Ugrid,const std::vector<double> Ve_grid,
	int sigma_type = 0){
	
	
	
	double mp = ME.at(element);
	double Q = QE.at(element);
	
	const auto &PM = BodyModel::Instance();
	if(!PM.isExist(element))
		return Function2<double>();

	MatrixElementType22 M22 = MET0(mp+mass);
	if(sigma_type == 1)
		M22 = MET1(mp+mass,mass*mp/(mp+mass)*U0);
	else if(sigma_type == -1)
		M22 = MET1(mp+mass,mass*mp/(mp+mass)*U0,1);
	else if(sigma_type == 2)
		M22 = MET2(mp+mass,mass*mp/(mp+mass)*U0,mass*mp/(mp+mass)*U0);
	
	return Function2<double>(Ugrid,Ve_grid,[Q,mp,mass,M22](double u,double ve){
			return sigmaC(mp,mass,sqrt(u*u+ve*ve),ve,M22,10,40)*Q*Q;
		});
	
	
}

#define KGeV 8.617333262e-14
Function2<double> SigmaElastic(double mass,const std::string &element,
	const std::vector<double> &Ugrid,const std::vector<double> &Ve_grid,
	int considerTemp = 0,
	int sigma_type = 0,size_t Nmk = 10000){
	
	
	const auto &PM = BodyModel::Instance();
	if(!PM.isExist(element))
		return Function2<double>();
	
	
	
	double mp = ME.at(element);
	
	Function1<double> WT(PM["Vesc"], apply_function<double,double>(
		PM["Temp"],[mp,considerTemp](double T){return considerTemp*sqrt(KGeV*T/mp);}));
	//std::cout << WT <<std::endl;
	
	return Function2<double>(Ugrid,Ve_grid,[WT,mp,mass,sigma_type,considerTemp,Nmk](double u,double ve){
			return sigmaTfacor(mp,mass,sqrt(u*u+ve*ve),ve,mp/(mp+mass)*U0,WT(ve),CAPTURE,sigma_type,Nmk);
		});
}

Function2<double> SigmaEscape(double mass,const std::string &element,
	const std::vector<double> v_grid,const std::vector<double> ve_grid,
	int sigma_type = 0,size_t Nmk = 10000){
	
	const auto &PM = BodyModel::Instance();
	if(!PM.isExist(element))
		return Function2<double>();
	
	
	double mp = ME.at(element);
	
	return Function2<double>(v_grid,ve_grid,[mp,mass,sigma_type,Nmk](double v,double ve){
			return sigmaTfacor(mp,mass,v,ve,mp/(mp+mass)*U0,0/*Temp*/,ESCAPE,sigma_type,Nmk);
		});
}

extern inline double thc(double x){
	double Ex = exp(x);
	if(x<0.5){
		return (1.0+x/2*(1.0+x/3*(1.0+x/4*(1.0+x/5*(1.+x/6*(1.+x/7*(1.+x/8*(1.+x/9*(1.+x/(10*(1.-x/11)))))))))))/Ex;
	}
	else
		return (Ex-1.0)/(Ex*x);
}

extern inline double fu(double u,double u0,double su){
	double D = 2*su*su;
	return exp(-(u-u0)*(u-u0)/(D))/pow(M_PI*D,1.5)*thc(4*u*u0/D);
}

double C_ND(Function1<double> Sigma,double su){
	auto U = Sigma.getXgrid();
	return integrateAB2(U,Sigma.getYgrid()*U*apply_function<double,double>(U,[su](double u){
			return 4*M_PI*U0*fu(u,U0,su);
		}) 
		);
}


std::vector<double> grid(size_t N,double Umin,double Umax,double q,double eps){
	std::vector<double> __grid(N);
	
	double et = pow(eps,1/q);
	double Z = pow(1+eps,1/q)-et; 
	
	double h = 1.0/(N-1);
	for(size_t i=0;i<N;i++){
		double t = i*h;
		__grid[i] = (pow(t*Z+et,q)-eps)*(Umax-Umin) + Umin;
	}
	return __grid;
}

std::string Filename(const std::string path,const std::string element,
				int isElastic,double mk,int sigmaType){
	std::string elastictype = "In";
	if(isElastic)
		elastictype = "El";
	return (path + elastictype + "(" + element +
		", " + std::to_string(mk) + ", " +  std::to_string(sigmaType) + ").dat");
}



#endif
