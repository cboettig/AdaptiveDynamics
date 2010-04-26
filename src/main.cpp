
int main(void)
{
	double sigma_mu = 0.05;
	double mu = 1.e-3;
	double sigma_c2 = .1; //gsl_pow_2(1.2);
	double sigma_k2 = 1;
	double ko = 1000;
	double xo = .5;

	int seed = time(NULL);

	double * phasetime;
	#pragma omp parallel private(phasetime) default(none) shared(sigma_mu, mu, sigma_k2, sigma_c2, ko, xo, seed)
	{
		phasetime = (double *) calloc(3,sizeof(double));
		branch_simulation(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, phasetime, &seed);
		free(phasetime);
	}


	int i;
	int npts = 50;
	double density[50];
	double times[50];
	for(i=0;i<npts;i++) times[i] = 1000.*i/(double) npts;

	analytics(&sigma_mu, &mu, &sigma_c2, &sigma_k2, &ko, &xo, times, density, &npts);

	return 0;
}

