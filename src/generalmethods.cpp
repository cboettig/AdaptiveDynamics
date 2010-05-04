#include "branch_vec.hh"
double sumrates(vector<pop> &poplist)
{
	vector<pop>::iterator p;
	double sum = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		sum += p->birth + p->death; 
	 }
	return sum;
}


double bdry(double x, par_list * pars){return -x*(1+pars->mc/pars->mk)/(1-pars->mc/pars->mk) ; }

int coexist(double x, double y, par_list * pars)
{ 
	if( (1 - C(x,y, pars)*K(x, pars)/K(y, pars) > EPSILON) && (1 - C(y,x, pars)*K(y, pars)/K(x, pars) > EPSILON) ) return 1;
	else return 0;
}

int branchcheck(vector<pop> &poplist, int threshold, par_list * pars)
{
	if( poplist.size() < 2 ){ return 0; }
	else{
		vector<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > threshold ){
				for(q=poplist.begin(); q != poplist.end(); q++){
					if(q->popsize > threshold && coexist(p->trait,q->trait, pars) ){
						return 1;
					}
				}
			}
		}
		return 0;
	}
}

int invade_pair(vector<pop> &poplist, int threshold, double * pair)
{
	if( poplist.size() < 3 ) return 0; 
	else {
		vector<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > threshold && (
				p->trait < GSL_MIN(pair[0], pair[1]) 
				|| p->trait > GSL_MAX(pair[0], pair[1]) 
				)){
//	printf("%lf, %lf, %lf %d\n",pair[0], pair[1], p->trait, p->popsize);
				return 1;
			}
		}
		return 0;
	}
}

int branches(vector<pop> &poplist, int threshold, double * pair, par_list * pars)
{
	if( poplist.size() < 2 ){ return 0; }
	else{
		vector<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > threshold ){
				for(q=poplist.begin(); q != poplist.end(); q++){
					if(q->popsize > threshold && coexist(p->trait,q->trait, pars) ){
						pair[0] = p->trait;
						pair[1] = q->trait;
						return 1;
					}
				}
			}
		}
		return 0;
	}
}
int mutual_invade(vector<pop> &poplist, int threshold, par_list * pars)
{
	vector<pop>::iterator p, q;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->popsize > threshold ){
			for(q=poplist.begin(); q != poplist.end(); q++){
				if(q->popsize > threshold && coexist(p->trait,q->trait, pars) ){
					double sq = (1 - C(p->trait,q->trait, pars)*K(p->trait, pars)/K(q->trait, pars) ), sp = (1 - C(p->trait,q->trait, pars)*K(q->trait, pars)/K(p->trait, pars) ) ;
					printf("%lf\t%lf\t", GSL_MIN(sp, sq), GSL_MAX(sp, sq));
//					if(p->popsize < q->popsize) p->parent_trait = 5; else q->parent_trait = 5; //flag the mutants
//					if(sp < sq) p->parent_trait = 5; else q->parent_trait = 5; //flag the weaker polymorphism
					return 1;
				}
			}
		}
	}
	return 0;
}


int finishline(vector<pop> &poplist, double line, int threshold)
{
	if( poplist.size() < 2 ){ return 0; }
	else{
		vector<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > threshold ){
				for(q=poplist.begin(); q != poplist.end(); q++){
					if(q->popsize > threshold && (p->trait-q->trait) > line ){
						return 1;
					}
				}
			}
		}
		return 0;
	}
}




int flagcheck(vector<pop> &poplist, double flag)
{
	vector<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->parent_trait == flag){return 1;}
	}
	return 0; 
}
double flagget(vector<pop> &poplist, double flag)
{
	vector<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->parent_trait == flag){return p->trait;}
	}
	return 0; 
}

void averagelist(vector<pop> &poplist, double sampletime, double *mean, int samplenumber)
{
	vector<pop>::iterator p;
	double meantrait = 0;
	double totalpop = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		totalpop += p->popsize;
		meantrait += p->popsize * p->trait;
	}
	mean[samplenumber] += meantrait/totalpop;
}

void printaverage(double *mean)
{
	int samplenumber;
	for(samplenumber=0; samplenumber < SAMPLES; samplenumber++){
		printf( "%lf, %lf\n", samplenumber*MAXTIME/SAMPLES, mean[samplenumber]/MAXTRIALS );
	}
}


void printlist(vector<pop> &poplist, double sampletime)
{
	vector<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		printf("%.2lf %d %lf %lf\n", sampletime, p->popsize, p->trait, p->parent_trait);
	}
	printf("\n");
}

void printfulllist(vector<pop> &poplist, double sampletime)
{
	vector<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		printf("t= %.2lf n= %d x= %lf p= %lf, b= %lf, d= %lf, c=  %lf\n", sampletime, p->popsize, p->trait, p->parent_trait, p->birth, p->death, p->compete);
	}
	printf("\n");
}

void printfreq(vector<pop> &poplist)
{
	vector<pop>::iterator p;
	p=poplist.begin();
	int N2 = p->popsize;
	double x2 = p->trait;
	if(poplist.size() == 2) p++;
	int N1  = p->popsize;
	double x1 = p->trait;
	printf("%d %8d %8lf %8d %8lf %8lf\n", N1, N2, (double) N1/( (double) N1+N2), N1+N2, x1, x2 );
}


