#include "branch_vec.hh"

int main(void)
{
	double sigma_mu = 0.03;
	double mu = 1.e-2;
	double sigma_c2 = .1; //gsl_pow_2(1.2);
	double sigma_k2 = 1;
	double ko = 100;
	double xo = .5;
	int threshold = 30;
	int seed = time(NULL);

	double * phasetime;
	#pragma omp parallel private(phasetime) default(none) shared(sigma_mu, mu, sigma_k2, sigma_c2, ko, xo, seed)
	{
		phasetime = (double *) calloc(6,sizeof(double));
		branch_simulation(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, phasetime, &seed, &threshold);
		free(phasetime);
	}

	double mean = 0;
	int i;
	int npts = 50;
	double density[50];
	double times[50];
	for(i=0;i<npts;i++) times[i] = 1000.*i/(double) npts;

//	analytics(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, times, density, &npts, &mean);
	printf("\n\n Mean: %g\n", mean );

	return 0;
}

