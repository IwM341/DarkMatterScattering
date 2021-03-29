#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "omp.h"

int main(void){
	
	omp_set_num_threads(1);
	int i;
	int N = 100000;
	clock_t t = clock();
	double *M = malloc(N*sizeof(double));
	
	int k = 0;
	
	while(k<100){
		for(int i = 0;i<N;i++){
			M[i] = (double)(rand())/RAND_MAX;
		}
		double sum =0;
		for(int i = 0;i<N;i++){
			sum += M[i];
		}
		printf("%lf\n",sum/N);
		k++;
	}

	printf("%d\n",(int)t);
	printf("%d\n",(int)clock());
	free(M);	
	return 0;
}