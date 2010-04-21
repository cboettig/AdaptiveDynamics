#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <list>

using namespace std;

typedef struct {
		int    popsize;
		double trait;
		double birth;
		double death;
		double parent_trait;
}pop ;


#define MK .5  // == 1/(2 sigma_k^2) ==> sigma_k^2 = 1
#define MC 2.5 // == 1/(2 sigma_c^2) ==> sigma_c^2 = .4
#define SIGMA_MU 0.05
#define MU 0.00005
#define Ko 1000
#define R 1.0
#define Xo .2
#define SAMPLES 1000 
#define MAXTIME 5e5
#define EPSILON 1e-4
#define MAXTRIALS 1

double C(double x, double y){return exp( -gsl_pow_2(x-y)*MC );}
double K(double x){return Ko*exp( -gsl_pow_2(x)*MK );}

double bdry(double x){return -x*(1+MC/MK)/(1-MC/MK) ; }
int coexist(double x, double y)
{ 
	if( 1 - C(x,y)*K(x)/K(y) > EPSILON && 1 - C(y,x)*K(y)/K(x) > EPSILON ) return 1;
	else return 0;
}
int branchcheck(list<pop> &poplist){
	list<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > 50 ){
				for(q=poplist.begin(); q != poplist.end(); q++){
					if(q->popsize > 50 && coexist(p->trait,q->trait) ){
						return 1;
					}
				}
			}
	}

	return 0;
}

void averagelist(list<pop> &poplist, double sampletime, double *mean, int samplenumber)
{
	list<pop>::iterator p;
	double meantrait = 0;
	double totalpop = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		totalpop += p->popsize;
		meantrait += p->popsize * p->trait;
	}
	mean[samplenumber] += meantrait/totalpop;
}

void printaverage(double *mean){
	int samplenumber;
	for(samplenumber=0; samplenumber < SAMPLES; samplenumber++){
		printf( "%lf, %lf\n", samplenumber*MAXTIME/SAMPLES, mean[samplenumber]/MAXTRIALS );
	}
}

void printlist(list<pop> &poplist, double sampletime)
{
	list<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		printf("%.2lf %d %lf %lf\n", sampletime, p->popsize, p->trait, p->parent_trait);
	}
	printf("\n");
}


void update_rates(list<pop> &poplist)
{
	list<pop>::iterator p, q;
	for(p=poplist.begin(); p != poplist.end(); p++){
		double compete = 0;
		for(q=poplist.begin(); q != poplist.end(); q++){
			compete += q->popsize * C( q->trait, p->trait );
		}

		p->birth = R* p->popsize;
		p->death = R* p->popsize * compete/K( p->trait );
	}

}


double sumrates(list<pop> &poplist)
{
	list<pop>::iterator p;
	double sum = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		sum += p->birth + p->death; 
	 }

	return sum;
}


void event(gsl_rng * rng, list<pop> &poplist, double sum)
{
	list<pop>::iterator p;
	double threshhold = sum*gsl_rng_uniform(rng);
	double cumulative = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		cumulative += p->birth;
		if( cumulative  > threshhold ){ 
			if(gsl_rng_uniform(rng) > MU){
				++p->popsize;
			} else {
				double step = gsl_ran_gaussian_ziggurat(rng,SIGMA_MU);
				pop pop_i = {1, p->trait+step, 0,0, p->trait};
				poplist.push_back(pop_i);
			}
			break;
		}
		cumulative += p->death;
		if( cumulative  > threshhold ){ 
			--p->popsize;
			if( p->popsize == 0 ){
				poplist.erase( p ); 
			}
			break;
		}
	}
}


int main(void)
{
	gsl_rng * rng = gsl_rng_alloc (gsl_rng_default); 
	gsl_rng_set(rng, time(NULL));

	list<pop> poplist;

	//double mean[SAMPLES]={0};

	printf(
		 "## MK = %lf, MC = %lf, SIGMA_MU = %lf, MU = %lf, Ko = %d, R = %lf, Xo = %lf\n",
		 MK, MC, SIGMA_MU, MU, Ko, R, Xo);

	int trial;
	for(trial = 0; trial < MAXTRIALS; trial++){

		double time = 0, sampletime = 0, Dt = MAXTIME/SAMPLES, sum;
		int samplenumber = 0;

		pop pop_i = {(int) round(K(Xo)), Xo, 0, 0, Xo};
		poplist.push_back(pop_i);

		int check = 1;
		while( time < MAXTIME ){
			update_rates(poplist);
			sum = sumrates(poplist);
			time += gsl_ran_exponential(rng, 1/sum);
			if(time > sampletime){
//				averagelist(poplist, sampletime, mean, samplenumber);
				printlist(poplist, sampletime);
				if(check){ if( branchcheck(poplist) ){ 
					check = 0;
					fprintf(stderr, "branched! %lf \n", time);
					//break;
				}}
				sampletime += Dt;
				++samplenumber;
			}
			event(rng, poplist, sum);
		}
//		printlist(poplist, sampletime);
		fprintf(stderr, "still branched = %d\n", branchcheck(poplist));		
		poplist.clear();
	}
	
//	printaverage(mean);
	gsl_rng_free(rng);
	return 0;
}
