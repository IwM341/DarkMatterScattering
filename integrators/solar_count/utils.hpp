#ifndef UTILS_H
#define UTILS_H
#include "csv_function.hpp"


#define U0 0.7667e-3



const std::map<string, double> ME({
			{"H1",1.0},{"He4",4.0},{"He3",3.0},{"C12",12.0},{"C13",13.0},
			{"N14",14.0},{"N15",15.0},{"O16",16.0},{"O17",17.0},{"O18",18.0},
			{"Ne",20.0},{"Na",23.0},{"Mg",24.0},{"Al",27.0},{"Si",28.0},{"P",31.0},
			{"S",32.0},{"Cl",35.5},{"Ar",40.0},{"K",39.0},{"Ca",40.0},
			{"Sc",45.0},{"Ti",48.0},{"V",51.0},{"Cr",52.0},{"Ca",40.0},
			{"Mn",55.0},{"Fe",55.8},{"Co",59.0},{"Ca",59.0}
			});
	
class SolarModel{
	std::map<std::string, std::vector<double>> SM;
	double Vesc;
	
	SolarModel():SM(Function::CSVTable("solar_model.dat")),
		Vesc(2.056e-3){
		
	}
	SolarModel(const SolarModel& root) = delete;
	SolarModel& operator = (const SolarModel&) = delete;
	
	
public:
	static SolarModel& Instance(){
		static SolarModel Inst();
		return Inst;
	}
	
	Function1<double> Row(const std::string & Title){
		if(Title == "Vesc"){
			auto X = SM["phi"];
			F.map([](double){return sqrt(x)*Vesc});
			return F;
		}
		return Function1(SM["Radius"],SM[column]);
	}
	
}

Function1<double> IntegrateSigma(Function2<double> SigmaVVesc,const std::string element){
	Function1<double> Sigma(SigmaVVesc.getXgrid());
	size_t N =  Sigma.size();
	
	auto SM = SolarModel::Instance();
	
	std::vector<double> VescS = SigmaVVesc.Ygrid(0);
	std::vector<double> US = SigmaVVesc.Xgrid();
	
	auto Vesc = Row("Vesc");
	auto Xi = Row("Radius");
	auto Rho = Row("RhoND");
	auto Element = Row(element);
	
	Function1<double> XiF(Vesc,Xi);
	Function1<double> RhoF(Vesc,Rho);
	Function1<double> ElementF(Vesc,Element);
	
	XiS = XiF(VescS);
	RhoS = RhoF(VescS);
	ElementS = ElementF(VescS);
	
	double mass = ME(element);
	
	#pragma omp parallel for 
	for(size_t i=0;i<N;i++){
		double U = US[i];
		Sigma.getY(i) = 3.0/(U0*U0*mass)*
		integrateAB2(VescS,
			apply_function(VescS,[vesc,U,&SigmaVVesc](double vesc){
				double v = sqrt(vesc*vesc+U*U)
				return v*v*SigmaVVesc(v,vesc);
				})*XiS*XiS*ElementS*RhoS
		);
	}
	return Sigma;
}



#endif
