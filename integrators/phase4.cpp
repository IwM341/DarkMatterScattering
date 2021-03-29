#include "phase4.hpp"

vec3::vec3():x(0),y(0),z(0){}

vec3::vec3(double x,double y,double z):
		x(x),y(y),z(z) {}

static vec3::vec3 Polar(double r,double theta,double phi){
	return vec3(r*cos(phi)*sin(theta),r*sin(phi)*sin(theta),r*cos(theta));
}


static vec3::vec3 PolarCos(double r,double cosTheta,double phi){
	double sinTheta = sqrt(1-cosTheta*cosTheta);
	return vec3(r*cos(phi)*sinTheta,r*sin(phi)*sinTheta,r*cosTheta);
}

vec3 vec3::operator +(const vec3 second) const{
	return vec3(x+second.x,y+second.y,z+second.z);
}

vec3 vec3::operator -(const vec3 second) const{
	return vec3(x-second.x,y-second.y,z-second.z);
}


vec3 vec3::operator -() const{
	return vec3(-x,-y,-z);
}

vec3 vec3::operator +=(const vec3 second){
	return vec3(x+=second.x,y+=second.y,z+=second.z);
}

vec3 vec3::operator -=(const vec3  second){
	return vec3(x-=second.x,y-=second.y,z-=second.z);
}

vec3 vec3::operator *(double a) const{
	return vec3(x*a,y*a,z*a);
}

vec3 vec3::operator /(double a) const{
	return vec3(x/a,y/a,z/a);
}

vec3 vec3::operator *=(double a){
	x*=a;
	y*=a;
	z*=a;
	return *this;
}

vec3 vec3::operator /=(double a){
	x/=a;
	y/=a;
	z/=a;
	return *this;
}

double vec3::operator*(vec3 second) const{
	return x*second.x+y*second.y+z*second.z;
}
