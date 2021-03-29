#ifndef PHASE_H
#define PHASE_H

#include <cmath>
#include <string>
class vec3{
public:
	double x,y,z;
	
	vec3();
	
	vec3(double x,double y,double z);
	
	static vec3 Polar(double r,double theta,double phi);
	
	static vec3 PolarCos(double r,double cosTheta,double phi);
	
	vec3 operator +(const vec3 second) const;
	vec3 operator -(const vec3 second) const;
	
	vec3 operator -() const;
	
	vec3 operator +=(const vec3 second);
	
	vec3 operator -=(const vec3 second);
	
	vec3 operator *(double a) const;
	vec3 operator /(double a) const;
	
	
	vec3 operator *=(double a);
	
	vec3 operator /=(double a);
	
	double operator*(vec3 second) const{
		return x*second.x+y*second.y+z*second.z;
	}
	
	double quad() const{
		return x*x+y*y+z*z;
	}
	
	double norm() const{
		return sqrt(quad());
	}
	
	operator std::string () const{
		return std::string("vec3(") + std::to_string(x) + "," + 
									std::to_string(y) + "," + 
									std::to_string(z) + ")";
	}
	

};

vec3 operator *(double a,vec3 P){
	return P*a; 
}

class vec4{
	public:
	double t;
	double x,y,z;
	
	vec4():t(0),x(0),y(0),z(0){}
	vec4(double t,double x,double y,double z):
		t(t),x(x),y(y),z(z) {}
	vec4(double t,vec3 X): t(t),x(X.x),y(X.y),z(X.z){}
	
	vec4(vec3 X,double mass):
		t(sqrt(mass*mass + X*X)),x(X.x),y(X.y),z(X.z){}
	
	vec3 vecPart() const{
		return vec3(x,y,z);
	}
	
	
	vec4 operator +(const vec4 second) const{
		return vec4(t+second.t,x+second.x,y+second.y,z+second.z);
	}
	vec4 operator -(const vec4 second) const{
		return vec4(t-second.t,x-second.x,y-second.y,z-second.z);
	}
	
	vec4 operator -() const{
		return vec4(-t,-x,-y,-z);
	}
	
	vec4 operator +=(const vec4 second){
		return vec4(t+=second.t,x+=second.x,y+=second.y,z+=second.z);
	}
	vec4 operator -=(const vec4 second){
		return vec4(t-=second.t,x-=second.x,y-=second.y,z-=second.z);
	}
	
	vec4 operator *(double a) const{
		return vec4(t*a,x*a,y*a,z*a);
	}
	vec4 operator /(double a) const{
		return vec4(t/a,x/a,y/a,z/a);
	}
	vec4 operator *=(double a){
		t*=a;
		x*=a;
		y*=a;
		z*=a;
		return *this;
	}
	vec4 operator /=(double a){
		t/=a;
		x/=a;
		y/=a;
		z/=a;
		return *this;
	}
	
	double operator*(vec4 second) const{
		return t*second.t - x*second.x-y*second.y-z*second.z;
	}
	
	double quad() const{
		return t*t- x*x - y*y - z*z;
	}
	
	operator std::string () const{
		return std::string("vec3(") + std::to_string(t) + "," +
									std::to_string(x) + "," + 
									std::to_string(y) + "," + 
									std::to_string(z) + ")";
	}
};

vec4 operator *(double a,vec4 P){
	return P*a; 
}

double E(double m,double p){
	return sqrt(m*m+p*p);
}
double P(double m,double E){
	return sqrt(E*E - m*m);
}
#endif
