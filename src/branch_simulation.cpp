#include "branch_vec.hh"
#include <omp.h>

/** Function for reporting phase of branching and stopping the simulation */
int checkphase(vector<pop> &poplist, char *PHASE, double pair[], double phasetime[], double sampletime, par_list * pars)
{	char phase = *PHASE;
	int done = 0;
	switch(phase)
	{
		case '0' : // first time two branches are established
			if(branches(poplist,pars->threshold, pair, pars) ){
				#if VERBOSE 
				printf("phase 1 done!\n"); printlist(poplist, sampletime);
				#endif
				phasetime[0] = sampletime;						
				phase = '2';
			}
			break;
			
		case '1' : // last time two branches restablished
			if(branches(poplist,pars->threshold, pair, pars) ){
				#if VERBOSE 
				printf("phase 1 done!\n"); printlist(poplist, sampletime);
				#endif
				phasetime[1] = sampletime;						
				phase = '2';
						}
			break;
			
		case '2' : // last time the dimorphic pair is invaded successfully
			if(!branchcheck(poplist,0, pars)){ 
				#if VERBOSE 
				printf("dimorphism lost! back to phase 1 !\n"); printlist(poplist, sampletime);
				#endif
				phase = '1'; 
				phasetime[4] += 1; // Number of times dimophism is lost in phase 2
			} else if(invade_pair(poplist, pars->threshold, pair)) {
				#if VERBOSE
				printf("phase 2 done!\n"); printlist(poplist, sampletime);
				#endif
				phasetime[2] = sampletime;
				phase = '3'; 
			}
			break;

		case '3' : 
			if(!branchcheck(poplist,0, pars)){ 
				#if VERBOSE
				printf("dimophism lost in phase 3!\n"); printlist(poplist, sampletime); 
				//printf("0 %.1lf %.1lf %.4lf %.4lf\n", phasetime[1] - phasetime[0], sampletime - phasetime[1], pair[0], pair[1]);
				#endif
				phase = '1'; 
				phasetime[5] += 1; // Number of times dimophism is lost in phase 1
			} else if( finishline(poplist, LINE, pars->threshold) ){
				#if VERBOSE
				printf("reached finishline!\n"); printlist(poplist, sampletime);
				// printf("1 %.1lf %.1lf %.4lf %.4lf\n", phasetime[1] - phasetime[0], sampletime - phasetime[1], pair[0], pair[1]); 
				#endif
				phasetime[3] = sampletime;
				phase = '4';
			}
			break;

		case '4' :
			done = 1;
			break;
	} // end switch
	*PHASE = phase;
	return done;
}



void initialize_population(vector<pop> &poplist, vector<CRow> &cmatrix, par_list * pars){

			double init_parent = M_PI; //initial id, arbitrary
			double N1o = round(K(pars->xo, pars));
			double comp = C(pars->xo,X2, pars)*N2o+N1o;
			double d = R*N1o*comp/K(pars->xo, pars);
			pop pop_i = {(int) N1o,	pars->xo, R*N1o, d, init_parent, comp};
			poplist.push_back(pop_i);
			grow_C(poplist, cmatrix, pars);
			if(N2o != 0){ /* Create the second clonal type if desired */
				comp = C(pars->xo,X2,pars)*N1o+N2o;
				d = R*N2o*comp/K(X2,pars);
				pop pop_j = {(int) N2o, X2, R*N2o,  d, init_parent, comp};
				poplist.push_back(pop_j);
				grow_C(poplist, cmatrix, pars);
			}
		}



/** Simulate a single replicate of evolutionary branching process.  Can record time to complete each phase of the process */
void branch_simulation(double *sigma_mu, double *mu, double *sigma_c2, double *sigma_k2, double *ko, double *xo, double * phasetime, int * seed, int *threshold)
{
	double mc = 1 / (2 * *sigma_c2);
	double mk = 1 / (2 * *sigma_k2);
	par_list p = {*sigma_mu, *mu, mc, mk, *ko, *xo, NULL, *threshold};
	par_list * pars = &p;

	gsl_rng *rng = gsl_rng_alloc (gsl_rng_default); 
	gsl_rng_set(rng, *seed);

	vector<pop> poplist;
	vector<CRow> cmatrix;

	
	
	/* Simulation */
	#pragma omp for private(poplist, cmatrix)
	for(int trial = 0; trial<MAXTRIALS; trial++){

		/* Create an initial population */
		initialize_population(poplist, cmatrix, pars);
		/* Reporting and tracking phase of branching */
		char phase = '0';
		double time = 0, sampletime = 0, Dt = MAXTIME/SAMPLES, sum;
		double pair[2]; 

		while( time < MAXTIME ){
			sum = sumrates(poplist);
			time += gsl_ran_exponential(rng, 1/sum);
			if(time > sampletime){
	//			printlist(poplist,sampletime);
				if(checkphase(poplist, &phase, pair, phasetime, sampletime, pars) ) break;
				sampletime += Dt;
			}
			event_and_rates(rng, poplist, sum, cmatrix, pars);
		}
		poplist.clear();
		cmatrix.clear();


		printf("%g %g %g %g %g %g\n", phasetime[0], phasetime[1], phasetime[2], phasetime[3], phasetime[4], phasetime[5]);
	}
	gsl_rng_free(rng);
}



void analytics(double *sigma_mu, double *mu, double *sigma_c2, double *sigma_k2, double *ko, double *xo, double *times, double *waiting_time_distribution, int * samples, double * mean)
{
	double mc = 1 / (2 * *sigma_c2);
	double mk = 1 / (2 * *sigma_k2);
	par_list p = {*sigma_mu, *mu, mc, mk, *ko, 1 / *ko, *xo, NULL};
	par_list * pars = &p;

	int i;
	for(i = 0; i < *samples; i++)
	{
		waiting_time_distribution[i] = waiting_time_1(times[i], pars);
//		printf("%g\n", waiting_time_distribution[i] );
	}
	*mean = mean_waiting_time_1(pars);

}

