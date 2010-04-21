#include "branch_vec.hh"
#include <omp.h>

int checkphase(vector<pop> &poplist, char *PHASE, double pair[], double phasetime[], double sampletime)
{	char phase = *PHASE;
	int done = 0;
	switch(phase)
	{
		case '1' :
			if(branches(poplist,THRESHOLD, pair) ){
				printf("phase 1 done!\n"); printlist(poplist, sampletime);
				phasetime[0] = sampletime;						
				phase = '2';
			}
			break;
			
		case '2' :
			if(!branchcheck(poplist,0)){ 
				printf("dimorphism lost! back to phase 1 !\n"); printlist(poplist, sampletime);
				phase = '1'; 
			} else if(invade_pair(poplist, THRESHOLD, pair)) {
				printf("phase 2 done!\n"); printlist(poplist, sampletime);
				phasetime[1] = sampletime;
				phase = '3'; 
			}
			break;

		case '3' : 
			if(!branchcheck(poplist,0)){ 
				printf("dimophism lost in phase 3!\n"); printlist(poplist, sampletime); 
				printf("0 %.1lf %.1lf %.4lf %.4lf\n", 
					phasetime[1] - phasetime[0], sampletime - phasetime[1], pair[0], pair[1]);
				phase = '1';
			} else if( finishline(poplist, LINE) ){
				printf("reached finishline!\n"); printlist(poplist, sampletime);
				printf("1 %.1lf %.1lf %.4lf %.4lf\n", 
					phasetime[1] - phasetime[0], sampletime - phasetime[1], pair[0], pair[1]); 
				printf("%lf\n", sampletime);
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

void branch_simulation(void)
{
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default); 
	gsl_rng_set(rng, time(NULL));

	vector<pop> poplist;
	vector<CRow> cmatrix;


	printf("##  MK = %.3lf, MC = %.3lf, SIGMA_MU = %.2e, MU = %.2e, Ko = %d, R = %.1lf \n", 
	MK, MC, SIGMA_MU, MU, Ko, R);
	printf(	"##  Xo = %.4lf, X2 = %.4lf N2o = %.0lf, Line = %.3lf continue = %d\n", Xo, X2, N2o, LINE, CONTINUE);
	printf("## Samples = %d, MaxTime = %.0e, MaxTrials = %d, Threshold = %d, epsilon = %.2e\n", 
	(int) SAMPLES, MAXTIME, MAXTRIALS, THRESHOLD, EPSILON);


	/* Initialize constants */
	double time = 0, sampletime = 0, Dt = MAXTIME/SAMPLES, sum;
	double init_parent = M_PI;

	/* Create an initial population */
	double N1o = round(K(Xo));
	double comp = C(Xo,X2)*N2o+N1o;
	double d = R*N1o*comp*Kinv(Xo);
	pop pop_i = {(int) N1o,	Xo, R*N1o, d, init_parent, comp};
	poplist.push_back(pop_i);
	grow_C(poplist, cmatrix);

	/* Reporting and tracking phase of branching */
	char phase = '1';
	double pair[2];
	double phasetime[2] = {0,0};

	while( time < MAXTIME ){
		sum = sumrates(poplist);
		time += gsl_ran_exponential(rng, 1/sum);

		if(time > sampletime){
			printlist(poplist,sampletime);
			
			if(checkphase(poplist, &phase, pair, phasetime, sampletime) ) break;

			sampletime += Dt;
		} // end sampling
		event_and_rates(rng, poplist, sum, cmatrix);
	} //end simulation
	poplist.clear();
	cmatrix.clear();
	gsl_rng_free(rng);
}

int main(void)
{
	branch_simulation();
	return 0;
}

