#include <cstdio>
#include "csv_function.hpp"
#include "integrator.hpp"

inline double earth_profile(const double rr)
// returns Earth's density profile in ton/m^3=g/cm^3
// r - in km
// from 1205.0930
// SD 08.07.2012
{
  const double r = rr*6371.0;
  double dens;
  if(r < 1217.1)
    {
      dens = 13.0122 - 8.45292*rr*rr;
    }
  else
    {
      if(r < 3458.7)
	{
	  dens = 12.5842 - 1.69929*rr - 1.94128*rr*rr
	    - 7.11215*rr*rr*rr;
	}
      else
	{
	  if(r < 5701.0)
	    {
	      dens = 6.8143 - 1.66273*rr - 1.18531*rr*rr;
	    }
	  else
	    {
	      if(r < 5951.0)
		{
		  dens = 11.1198 - 7.87054*rr;
		}
	      else
		{
		  if(r < 6336.0)
		    {
		      dens = 7.15855 - 3.85999*rr;
		    }
		  else
		    {
		      if(r < 6351.0)
			{
			  dens = 2.92;
			}
		      else // rr < 6371.0
			{
			  dens = 2.72;
			}
		    }
		}
	    }
	}
    }
  return dens; // in g/cm^3 = ton/m^3

}
int main(void){
	double rho_av = integrateAB5([](double ksi){
			return earth_profile(ksi)*3*ksi*ksi;
		},0,1,1000);
	std::function<double(double)> rho = [rho_av](double x){return earth_profile(x)/rho_av;};
	
	std::function<double(double)> mass = [rho](double m){return integrateAB5([rho](double x){
			return rho(x)*3*x*x;
		},0,m,1000);};
	
	
	std::vector<double> Xi(201);
	double h = 1.0/(Xi.size()-1);
	for(size_t i=0;i<Xi.size();i++)
		Xi[i] = h*i;
	
	Function1<double> Rho(Xi,apply_function<double,double>(Xi,rho));
	Function1<double> M(Xi,apply_function<double,double>(Xi,mass));
	std::ofstream out("earth.dat");
	out << Rho;
	std::ofstream out1("earthM.dat");
	out1 << M;
	
}
