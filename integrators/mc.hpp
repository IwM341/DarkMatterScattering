#ifndef MC_INT
#define MC_INT

#include <cstdlib>
#include "omp.h"
#include <functional>
#include <utility>

template<typename CoordType>
class AreaGenerator{
	public:
	virtual CoordType RandomCoords() const = 0 ;
	virtual double AreaVolume() const = 0;
};

class MonteCarloIntegrator{
	public:	
	template<typename ReturnFunctionType,typename CoordType>
	ReturnFunctionType operator() (
			std::function<ReturnFunctionType(CoordType)> function,
			AreaGenerator<CoordType> *Gen,int N = 2)
	{
		
		ReturnFunctionType sum = 0;
		#pragma omp parallel for reduction(+:sum) shared(N)
		for(int i=0;i<N;i++){
			sum += function(Gen->RandomCoords());
		}
		return sum*Gen->AreaVolume()/N;
	}
};

#endif
