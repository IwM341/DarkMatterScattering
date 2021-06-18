#ifndef UTILS_H
#define UTILS_H
#include "csv_function.hpp"
#include "integrator.hpp"
#include "cross_section3.hpp"

#define U0 0.7667e-3



const std::map<std::string, double> ME({
			{"H1",1.0},{"He4",4.0},{"He3",3.0},{"C12",12.0},{"C13",13.0},
			{"N14",14.0},{"N15",15.0},{"O16",16.0},{"O17",17.0},{"O18",18.0},
			{"Ne",20.0},{"Na",23.0},{"Mg",24.0},{"Al",27.0},{"Si",28.0},{"P",31.0},
			{"S",32.0},{"Cl",35.5},{"Ar",40.0},{"K",39.0},{"Ca",40.0},
			{"Sc",45.0},{"Ti",48.0},{"V",51.0},{"Cr",52.0},{"Ca",40.0},
			{"Mn",55.0},{"Fe",55.8},{"Co",59.0},{"Ca",59.0}
			});
	
class BodyModel{
	std::map<std::string, std::vector<double>> BM;
	double Vesc;
	
	BodyModel():BM(Function::CSVTable("earth_model.dat")),
		Vesc(3.7336e-5){
		
	}
	BodyModel(const BodyModel& root);
	BodyModel& operator = (const BodyModel&);
	
	
public:
	static BodyModel &Instance(){
		static BodyModel Inst;
		return Inst;
	}
	
	std::vector<double> operator[](const std::string & column) const{
		if(column == "Vesc"){
			double vesc = Vesc;
			return apply_function<double,double>(BM.at("phi"),[vesc](double x)->double{return sqrt(x)*vesc;});
		}
		return BM.at(column);
	}
	double v_esc() const{
		return Vesc;
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
		Sigma.getY(i) = 3.0/(U0*U0*mass)*
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
	
	double mp = ME[element];
	
	const auto &PM = BodyModel::Instance();
	
	MatrixElementType23 M23 = MET0Q(mp+mass);
	if(sigma_type == 1)
		M23 = MET1Q(mp+mass,mass*mp/(mp+mass)*U0);
	else if(sigma_type == 2)
		M23 = MET1Q(mp+mass,mass*mp/(mp+mass)*U0,1);
	else if((sigma_type == 3)
		M23 = MET2Q(mp+mass,mass*mp/(mp+mass)*U0,mass*mp/(mp+mass)*U0);
	
	return Function2<double> SgIn(Ugrid,Ve_grid,[](double u,double ve){
			return sigmaC(ME[element] mass,sqrt(u*u+ve*ve),ve,M23,2,2,4,40);
		});
	
	
}

#endif
