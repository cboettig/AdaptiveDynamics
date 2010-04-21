#include "branch_vec.hh"
#include <omp.h>


void branch_simulation(void)
{
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default); 
	gsl_rng_set(rng, time(NULL));
	int trial;

	vector<pop> poplist;
	vector<CRow> cmatrix;
	double mean[SAMPLES]={MAXTIME};



	printf("##  MK = %.3lf, MC = %.3lf, SIGMA_MU = %.2e, MU = %.2e, Ko = %d, R = %.1lf \n", 
	MK, MC, SIGMA_MU, MU, Ko, R);
	printf(	"##  Xo = %.4lf, X2 = %.4lf N2o = %.0lf, Line = %.3lf continue = %d\n", Xo, X2, N2o, LINE, CONTINUE);
	printf("## Samples = %d, MaxTime = %.0e, MaxTrials = %d, Threshold = %d, epsilon = %.2e\n", 
	(int) SAMPLES, MAXTIME, MAXTRIALS, THRESHOLD, EPSILON);


/*
    // * Section for looping over different parameters * /
 	int i, j;
	double Xo;
	#pragma omp parallel for default(none) shared(rng) private(poplist, cmatrix, Xo,i,j,trial)
	for(i = 0; i < 40; i++){
		Xo = -.4 + i*0.02;
		double X2 = -.2;
		for(j=0; j<40; j++){
			X2 = -.4 + j*0.02;
			double N2o = round(K(X2));
*/


		#pragma omp parallel for default(none) shared(rng, mean) private(poplist, cmatrix, trial)
		for(trial = 0; trial < MAXTRIALS; trial++){

			/* generating png files */
			char filename[255];
			sprintf(filename, 
			    "plots/trial%d_N%d_MK%.2lf_MC%.0lf_SMU%.1e_MU%.1e.png",
				trial, Ko, MK, MC, SIGMA_MU, MU);
			#if PNG_ON
			pngwriter png(X_PLOTSIZE, Y_PLOTSIZE, 1.0, filename);
			#endif


			/* Initialize constants */
			double time = 0, sampletime = 0, Dt = MAXTIME/SAMPLES, sum;
			int samplenumber = 0;
			double init_parent = M_PI;

			/* Create an initial population */
			double N1o = round(K(Xo));
			double comp = C(Xo,X2)*N2o+N1o;
			double d = R*N1o*comp*Kinv(Xo);
			pop pop_i = {(int) N1o,	Xo, R*N1o, d, init_parent, comp};
			poplist.push_back(pop_i);
			grow_C(poplist, cmatrix);

			if(N2o != 0){ /* Create the second starting clonal type if desired */
				comp = C(Xo,X2)*N1o+N2o;
				d = R*N2o*comp*Kinv(X2);
				pop pop_j = {(int) N2o, X2, R*N2o,  d, init_parent, comp};
				poplist.push_back(pop_j);
				grow_C(poplist, cmatrix);
			}

			/* Reporting and tracking phase of branching */
			double phase1_time = 0, phase2_time=0;
			double pair[2];
			// To get coexist stats, we want both of these off.  For branching times, check should start at 1!!!
			int phase1 = 1, phase2 = 0, phase3 = 0;


			while( time < MAXTIME ){
				sum = sumrates(poplist);
				time += gsl_ran_exponential(rng, 1/sum);
/*
				if( poplist.size() < 2 ){
					printf("%.1lf %.2lf, %.2lf, %d\n", time, Xo, X2, poplist.begin()->trait == Xo); 
					break;
				}
*/
				if(time > sampletime){
//					averagelist(poplist, sampletime, mean, samplenumber);
//					printlist(poplist,sampletime);	
					#if PNG_ON					
					plotpng(poplist, (int) time/Dt, png);
					#endif	

					if(phase1){if(branches(poplist,THRESHOLD, pair) ){
						if(VERBOSE){printf("phase 1 done!\n"); printlist(poplist, sampletime);}
						phase1_time = sampletime;						
						phase2 = 1; phase1 = 0;
					}}
					if(phase2){if(!branchcheck(poplist,0)){ 
							if(VERBOSE){printf("dimorphism lost! back to phase 1 !\n"); printlist(poplist, sampletime);}
							phase1= 1; phase2=0;					
					}}

					if(phase2){if(invade_pair(poplist, THRESHOLD, pair)){
						if(VERBOSE){printf("phase 2 done!\n"); printlist(poplist, sampletime);}
						phase2_time = sampletime;
						phase3 = 1; phase2 = 0;
//						dimorph = gettraits(poplist);
//						int n_types = poplist.size();
					}}
					if(phase3){if(!branchcheck(poplist,0)){ 
							if(VERBOSE){ printf("dimophism lost in phase 3!\n"); printlist(poplist, sampletime); 
								printf("0 %.1lf %.1lf %.4lf %.4lf\n", 
									phase2_time - phase1_time, sampletime - phase2_time, pair[0], pair[1]);
							}
								phase3 = 0; phase1 = 1;
							if(!CONTINUE) break;					
					}}
			         if(phase3){if( finishline(poplist, LINE) ){
	                        if(VERBOSE) { printf("reached finishline!\n"); printlist(poplist, sampletime);
								printf("1 %.1lf %.1lf %.4lf %.4lf\n", 
									phase2_time - phase1_time, sampletime - phase2_time, pair[0], pair[1]); 
							}
							printf("%lf\n", sampletime);
	                        break;
                    }}

					sampletime += Dt;
					++samplenumber;
				}

				event_and_rates(rng, poplist, sum, cmatrix);
			}
			mean[trial] = sampletime;
			poplist.clear();
			cmatrix.clear();
			#if PNG_ON
			png.close();
			#endif
		}
		printf("%lf, %lf\n", gsl_stats_mean(mean, 1, MAXTRIALS), sqrt( gsl_stats_variance(mean, 1, MAXTRIALS) ));
//	}}
	gsl_rng_free(rng);
}

int main(void)
{
	branch_simulation();
	return 0;
}

