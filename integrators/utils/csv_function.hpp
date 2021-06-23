#ifndef CSV_FUNCTION_H
#define CSV_FUNCTION_H

#include <functional>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility> 
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <map>

template <class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& V){
	os << "Vector[ " ;
	for (size_t i=0;i<V.size();i++) {
		if(i < V.size()-1)
			os << V[i] << ", ";
		else
			os << V[i] << "]";
    }
	return os;
}

template <typename T1,typename T2>
std::vector<T2> apply_function(const std::vector<T1> &X,std::function<T2(T1)> F){
	std::vector<T2> Y(X.size());
	for(size_t i=0;i<X.size();i++){
		Y[i] = F(X[i]);
	}
	return Y;
}

template <typename T>
std::vector<T> operator +(const std::vector<T> & X,const std::vector<T> & Y){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]+Y[i];
	return Z;
}


template <typename T>
std::vector<T> operator *(const std::vector<T> & X,const std::vector<T> & Y){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]*Y[i];
	return Z;
}

template <typename T>
std::vector<T> operator *(const std::vector<T> & X,T Y){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]*Y;
	return Z;
}

template <typename T>
std::vector<T> operator *(T Y,const std::vector<T> & X){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]*Y;
	return Z;
}

template <typename T>
std::vector<T> operator -(const std::vector<T> & X,const std::vector<T> & Y){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]-Y[i];
	return Z;
}

template <typename T>
std::vector<T> operator /(const std::vector<T> & X,const std::vector<T> & Y){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]/Y[i];
	return Z;
}

template <typename T>
std::vector<T> operator /(const std::vector<T> & X,T Y){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = X[i]/Y;
	return Z;
}

template <typename T>
std::vector<T> operator /(T Y,const std::vector<T> & X){
	std::vector<T> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = Y/X[i];
	return Z;
}


extern std::vector<double> sqrt(const std::vector<double> &X){
	std::vector<double> Z(X.size());
	for(size_t i=0;i<X.size();i++)
		Z[i] = sqrt(X[i]);
	return Z; 
}

template <typename T>
std::vector<T> parse_string(const std::string &S){
	std::stringstream ss(S);
	std::vector<T> parsed;
	T it;
	while(ss>>it)
		parsed.push_back(it);
	return parsed;
}



extern inline find_less(const std::vector<double> &X,double x,
						size_t i1 = 0,size_t i2 = -1){
	size_t N = X.size();
	if(i2 == -1)
		i2 = N-1; 
	
	while(i1 < i2-1){
		size_t i = (i1 + i2)/2;
		if(x < X[i]){
			i2 = i;
		}
		else if(x > X[i]){
			if(x < X[i1+1])
				return i1;
			i1 = i;
		}
		else
			return i;
	}
	return i1;
}


template <class T>
class Function1{
	//std::vector<std::pair<double,T>> XY;
	std::vector<double> X;
	std::vector<double> Y;
	friend std::ostream& operator<<(std::ostream& os, const Function1<T>& F){
		
		size_t prec = os.precision();
		os <<std::setprecision(17);
		for (size_t i=0;i<F.size();i++) {
			os << F.X[i] << '\t' <<  F.Y[i]<< std::endl;
		}
		os << std::setprecision(prec);
		return os;
	}
	void sort(){
		std::vector<std::pair<double,T>> XY(X.size());
		for(int i=0;i<X.size();i++){
			XY[i] = std::pair<double,T>(X[i],Y[i]);
		}
		std::sort(XY.begin(),XY.end(),
			[](std::pair<double,T> a,std::pair<double,T> b){
				return a.first<b.first;
			}
		);
		for(int i=0;i<X.size();i++){
			X[i] = XY[i].first;
			Y[i] = XY[i].second;
		}
	}
	friend Function1<double> Swap(const Function1<double> &F){
		return Function1<double>(F.Y,F.X);
	}
public:
	Function1(size_t N = 0):X(N),Y(N){}
	Function1(const std::vector<double> &X):X(X),Y(X.size()){}
	Function1(const std::vector<double> &X,std::function<T(double)> f,bool sorted = false):
		X(X){
			if(!sorted)
				std::sort(this->X.begin(),this->X.end());
			Y = apply_function<double,T>(this->X,f);
		}
	Function1(std::vector<std::pair<double,T>> XY,bool sorted = false){
		for(auto xy : XY){
			X.push_back(xy.first);
			Y.push_back(xy.second);
		}
		if(!sorted)
			sort();
	}
	Function1(const std::vector<double> &X,const std::vector<T> &Y,bool sorted = false):
	X(X),Y(Y){if(!sorted) sort();}
	

	double &getX(size_t i){return X[i];}
	double &getY(size_t i){return Y[i];}
	
	std::vector<double> getXgrid() const{
		return X;
	}
	std::vector<double> getYgrid() const{
		return Y;
	}
	
	void push_back(double x,double y){
		X.push_back(x);
		Y.push_back(y);
		sort();
	}
	void push_back(const std::vector<double> &x,const std::vector<double> &y){
		if(x.size()!=y.size()){
			std::cout << "error: x.size() != y.size()" << std::endl;
		}
		else{
			X.insert(X.end(),x.begin(),x.end());
			Y.insert(Y.end(),y.begin(),y.end());
		}
		sort();
	}
	void insert(Function1<T> F){
		push_back(F.x,F.Y);
	}
	T UniformF(double x){
		size_t i = (x-X[0])/(X[1]-X[0]);
		size_t i1 = i+1;
		if(i1 >= X.size()){
			return Y[X.size()-1];
		}
			
		double a = (x-X[i])/(X[i1]-X[i]);
		return (1-a)*Y[i]+a*Y[i1];
	}
	T operator ()(double x) const{
		size_t N = X.size();
		size_t i = find_less(X,x);
		size_t i1 = i+1;
		if(i1 >= N){
			return Y[i];
		}
		double a = (x-X[i])/(X[i1]-X[i]);
		return (1-a)*Y[i]+a*Y[i1];
	}
	
	void map(std::function<T(T)> f){
		for(auto & y : Y){
			y = f(y);
		}
	}
	
	size_t size() const{return X.size();}

};

namespace Function{
	Function1<double> LoadFunction1(const char *filename){
		Function1<double> F;
		std::ifstream ifs(filename, std::ifstream::in);
		double x;
		double y;
		while(ifs>>x>>y){
			F.push_back(x,y);
		}
		return F;
	}
	
}

template <class T>
class Function2{
	std::vector<double> X;
	std::vector<double> Y;
	std::vector<T> F;
	std::vector<size_t> ny;
	std::vector<size_t> Ny;
	
	friend std::ostream& operator<<(std::ostream& os, const Function2<T>& F){
		
		size_t prec = os.precision();
		os <<std::setprecision(17);
		for (size_t i=0;i<F.X.size();i++) {
			os << F.X[i];
			for(size_t j = 0;j<F.ny[i];j++){
				os <<'\t' << F.getY(i,j) << '\t' << F.getF(i,j) << std::endl;
			}
		}
		os << std::setprecision(prec);
		return os;
	}
	void fullNy(){
		Ny.resize(ny.size()+1);
		Ny[0] = 0;
		for(size_t i=0;i<ny.size();i++){
			Ny[i+1] = Ny[i]+ny[i];
		}
	}
public:
	

	void fullF(std::function<T(double,double)> f){
		if(F.size() !=Y.size())
			F.resize(Y.size());
		#pragma omp parallel for 
		for(size_t i = 0;i<X.size();i++){
			for(size_t j = 0;j<ny[i];j++){
				F[Ny[i] + j] = f(X[i],Y[Ny[i] + j]);
			}
		}
	}
	Function2(){}
	
	Function2(std::vector<double> X,std::vector<double> Y,
				std::vector<size_t> ny,std::vector<T> F):X(X),Y(Y),F(F),ny(ny){
		fullNy();
	}
	Function2(std::vector<double> X,std::vector<double> Y,
				std::vector<size_t> ny,std::function<T(double,double)> f):
				X(X),Y(Y),ny(ny),F(Y.size()){
		fullNy();
		fullF(f);
		
	}
	Function2(std::vector<double> Xgrid,
				std::vector<double> Ygrid,std::vector<T> F ):X(Xgrid),
				Y(Xgrid.size()*Ygrid.size()),F(F),
				ny(Xgrid.size(),Ygrid.size()){
		fullNy();
		size_t _Ny = Ygrid.size();
		for(size_t i=0;i<ny.size();i++){
			for(size_t j=0;j<Ygrid.size();j++)
				Y[_Ny*i+j] = Ygrid[j];
		}
	}
	Function2(std::vector<double> Xgrid,
				std::vector<double>Ygrid, std::function<T(double,double)> f):X(Xgrid),
				Y(Xgrid.size()*Ygrid.size()),F(Xgrid.size()*Ygrid.size()),
				ny(Xgrid.size(),Ygrid.size()){
		fullNy();
		
		size_t _Ny = Ygrid.size();
		for(size_t i=0;i<ny.size();i++){
			for(size_t j=0;j<Ygrid.size();j++)
				Y[_Ny*i+j] = Ygrid[j];
		}
		fullF(f);
	}
	
	Function2(std::vector<double> Xgrid,
				std::vector<double>Ygrid):X(Xgrid),
				Y(Xgrid.size()*Ygrid.size()),F(Xgrid.size()*Ygrid.size()),
				ny(Xgrid.size(),Ygrid.size()){
		fullNy();
		
		size_t _Ny = Ygrid.size();
		for(size_t i=0;i<ny.size();i++){
			for(size_t j=0;j<Ygrid.size();j++)
				Y[_Ny*i+j] = Ygrid[j];
		}
	}
	
	Function2(size_t Nx,size_t Ny):X(Nx),Y(Nx*Ny),F(Nx*Ny),
		ny(Nx,Ny),Ny(Nx+1){
		fullNy();
	}
	Function2(size_t Nx,std::vector<size_t> ny):X(Nx),
		ny(ny),Ny(Nx+1){
		fullNy();
		Y.resize(Ny[Nx]);
		F.resize(Y.size());
	}
	
	T operator () (double x,double y) const{
		size_t i1 = find_less(X,x);
		std::cout << "i1 = " <<i1 <<std::endl; 
		if(i1 >= X.size()-1){
			
			size_t j1 = find_less(Y,y,Ny[i1]);
			if(j1 == Y.size())
				return F[j1];
			else{
				double a = (y-Y[j1])/(Y[j1+1]-Y[j1]);
				return (1-a)*F[j1] + a*F[j1+1];
			}
		}
		else{
			size_t j1 = find_less(Y,y,Ny[i1],Ny[i1+1]-1);
			size_t j2 = find_less(Y,y,Ny[i1+1],Ny[i1+2]-1);
			//std::cout << i1 << " " << j1 << " " <<j2 << std::endl;
			double f0y;
			double fhy;
			if(j1 >= Ny[i1+1]-1)
				f0y = F[j1];
			else{
				double ay0 = (y-Y[j1])/(Y[j1+1]-Y[j1]);
				f0y = (1-ay0)*F[j1] + ay0*F[j1+1];
			}
			//std::cout << f0y <<std::endl;
			if(j2 >= Ny[i1+2]-1)
				fhy = F[j1];
			else{
				double ayh = (y-Y[j2])/(Y[j2+1]-Y[j2]);
				fhy = (1-ayh)*F[j2] + ayh*F[j2+1];
			}
			//std::cout << fhy <<std::endl;
			double axh = (x-X[i1])/(X[i1+1]-X[i1]);
			return (1-axh)*f0y + axh*fhy;
		}
		
	}
	std::vector<double> getXgrid()const{return X;}
	std::vector<double> getYgrid(size_t i) const{
			return std::vector<double> (Y.begin() + Ny[i],Y.begin() + Ny[i+1]);
		}
	double &getX(size_t i){
		return X[i];
	}
	double &getY(size_t i,size_t j){
		return Y[Ny[i] + j];
	}
	T &getF(size_t i,size_t j){
		return F[Ny[i] + j];
	}
	const double getX(size_t i) const {
		return X[i];
	}
	const double getY(size_t i,size_t j) const {
		return Y[Ny[i] + j];
	}
	const T getF(size_t i,size_t j) const{
		return F[Ny[i] + j];
	}
	
	void map(std::function<T(double)> f){
		std::for_each(F.begin(), F.end(), f);
	}
	
	std::string toString() const{
		std::stringstream ss;
		ss << X << std::endl;
		ss << ny << std::endl;
		ss << Ny << std::endl;
		ss << Y << std::endl;
		ss << F << std::endl;
		
		return ss.str(); 
	}
	
	Function1<T> atX(size_t i) const{
		return Function1<T>(std::vector<double>(Y.begin()+Ny[i],Y.begin()+Ny[i+1]),
						std::vector<T>(F.begin()+Ny[i],F.begin()+Ny[i+1]),true);
	}
};

namespace Function{
	Function2<double> LoadFunction2(const char *filename){
		std::vector<double> X;
		std::vector<double> Y;
		std::vector<double> F;
		std::vector<size_t> ny;
		std::vector<size_t> Ny;
		std::ifstream ifs(filename, std::ifstream::in);
		
		std::vector<double> nums;
		while(!ifs.eof()){
			std::string S;
			std::getline(ifs,S);
			nums = parse_string<double>(S);
			if(nums.size() == 3){
				X.push_back(nums[0]);
				Y.push_back(nums[1]);
				F.push_back(nums[2]);
				ny.push_back(1);		
			}
			else if(nums.size() == 2){
				Y.push_back(nums[0]);
				F.push_back(nums[1]);
				ny.back() += 1;
			}
		}
		return Function2<double>(X,Y,ny,F);
	}
}


namespace Function{
	
extern inline std::map<std::string, std::vector<double>> CSVTable(const char * filename){
	std::map<std::string, std::vector<double>> Funcs;
	std::ifstream ifs(filename, std::ifstream::in);
	
	if(!ifs.is_open()){
		std::cout << "CSVTable error: no such file" <<std::endl;
		return std::map<std::string, std::vector<double>>();
	}
	
	std::string S;
	std::getline(ifs,S);
	std::vector<std::string> cols = parse_string<std::string>(S);
	
	for(auto title : cols){
		Funcs[title] = std::vector<double>();
	}
	
	std::vector<double> nums;
	while(!ifs.eof()){
		std::getline(ifs,S);
		nums = parse_string<double>(S);
		if(nums.size()>=cols.size()){
			for(size_t i=0;i<cols.size();i++){
				Funcs[cols[i]].push_back(nums[i]);
			}
		}
	}
	return Funcs;
}

}
#endif
