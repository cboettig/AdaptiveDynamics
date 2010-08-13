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

	double * phasetime, * xpair, * ypair;
	#pragma omp parallel private(phasetime) default(none) shared(sigma_mu, mu, sigma_k2, sigma_c2, ko, xo, seed)
	{
		phasetime = (double *) calloc(6,sizeof(double));
		xpair = (double *) calloc(6,sizeof(double));
		ypair = (double *) calloc(6,sizeof(double));
		branch_simulation(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, phasetime, &seed, &threshold, xpair, ypair);
		printf("\n");
		printf(" %g time to first dimorphism\n %g time persistant dimorphism was established\n %g time dimorphism invaded\n %g time at finish pt\n %g collapses after dimorphism invaded\n %g collapses after dimorphism established\n\n", phasetime[0], phasetime[1], phasetime[2], phasetime[3], phasetime[4], phasetime[5]);
		int i;
		printf("xpair \t\t ypair, \t\t diff\n");
		for(i=0;i<4;i++) printf("%f\t %f\t %f\n", xpair[i], ypair[i], xpair[i]-ypair[i]);		
		free(phasetime);
		free(xpair);
		free(ypair);
	}

	double mean = 0;
	int i;
	int npts = 50;
	double density[50];
	double times[50];
	for(i=0;i<npts;i++) times[i] = 1000.*i/(double) npts;

	/* defined in branch_simulation w/ fns from ecoevo_model.cpp */
	analytics(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, times, density, &npts, &mean);
	printf("\n\n Analytic Mean: %g\n", mean );


/* defined in analytics.cpp, rather slow.  quality of coexist approx tbd */
//	analytic_contours_wrapper(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo);


	return 0;
}

