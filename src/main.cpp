#include "branch_vec.hh"

int main(void)
{
	double maxtime = 5e4;
	int samples = 1e4;
	double sigma_mu = 0.03;
	double mu = 1.e-2;
	double sigma_c2 = .1; //gsl_pow_2(1.2);
	double sigma_k2 = 1;
	double ko = 100;
	double xo = .1;
	int threshold = 30;
	int seed = time(NULL);

	double * phasetime, * xpair, * ypair;


	printf("Testing branch_simulation...\n");
	#pragma omp parallel private(phasetime, xpair, ypair) default(none) shared(sigma_mu, mu, sigma_k2, sigma_c2, ko, xo, seed, threshold, samples, maxtime)
	{
		phasetime = (double *) calloc(6,sizeof(double));
		xpair = (double *) calloc(6,sizeof(double));
		ypair = (double *) calloc(6,sizeof(double));
		branch_simulation(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, phasetime, &seed, &threshold, xpair, ypair, &maxtime, &samples);
		printf("\n");
		printf(" %g time to first dimorphism\n %g time persistant dimorphism was established\n %g time dimorphism invaded\n %g time at finish pt\n %g collapses after dimorphism invaded\n %g collapses after dimorphism established\n\n", phasetime[0], phasetime[1], phasetime[2], phasetime[3], phasetime[4], phasetime[5]);
		int i;
		printf("xpair \t\t ypair, \t\t diff\n");
		for(i=0;i<4;i++) printf("%f\t %f\t %f\n", xpair[i], ypair[i], xpair[i]-ypair[i]);		
		free(phasetime);
		free(xpair);
		free(ypair);
	}

	int i, npts = 50;
	double mean=0, density[50], times[50];
	
	for(i=0;i<npts;i++) times[i] = 1000.*i/(double) npts;

	/* defined in branch_simulation w/ fns from ecoevo_model.cpp */
	printf("Testing analytic mean time to branching approximation...\n");
	analytics(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, times, density, &npts, &mean);
	printf("\n\n Analytic Mean: %g\n", mean );

	double coexist_time[400], xval[400], yval[400];


	/* coexistence time */
	printf("Testing coexist_simulation...\n");
	coexist_simulation(	&sigma_mu, &mu, &sigma_c2, &sigma_k2, 
						&ko, &xo, coexist_time, &seed, &threshold, 
						xval, yval, &maxtime, &samples);
	for(i=0;i<400;i++) printf("%g %g %g\n", coexist_time[i], xval[i], yval[i]);

/* defined in analytics.cpp, rather slow.  quality of coexist approx tbd */
printf("Testing analytic_contours_wrapper...\n");
	analytic_contours_wrapper(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, coexist_time, xval, yval);
	for(i=0;i<400;i++) printf("%g %g %g\n", coexist_time[i], xval[i], yval[i]);

printf("Done!\n");
	return 0;
}

