#include "branch_vec.hh"

double C(double x, double y){return exp( -gsl_pow_2(x-y)*MC );}
double K(double x){return Ko*exp( -gsl_pow_2(x)*MK );}
double Kinv(double x){return IK*exp( gsl_pow_2(x)*MK );}

double sumrates(vector<pop> &poplist)
{
	vector<pop>::iterator p;
	double sum = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		sum += p->birth + p->death; 
	 }

	return sum;
}


double bdry(double x){return -x*(1+MC/MK)/(1-MC/MK) ; }
int coexist(double x, double y)
{ 
	if( (1 - C(x,y)*K(x)/K(y) > EPSILON) && (1 - C(y,x)*K(y)/K(x) > EPSILON) ) return 1;
	else return 0;
}

int branchcheck(vector<pop> &poplist, int threshold)
{
	if( poplist.size() < 2 ){ return 0; }
	else{
		vector<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > threshold ){
				for(q=poplist.begin(); q != poplist.end(); q++){
					if(q->popsize > threshold && coexist(p->trait,q->trait) ){
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

int branches(vector<pop> &poplist, int threshold, double * pair)
{
	if( poplist.size() < 2 ){ return 0; }
	else{
		vector<pop>::iterator p, q;
		for(p=poplist.begin(); p != poplist.end(); p++){
			if(p->popsize > threshold ){
				for(q=poplist.begin(); q != poplist.end(); q++){
					if(q->popsize > threshold && coexist(p->trait,q->trait) ){
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
int mutual_invade(vector<pop> &poplist, int threshold)
{
	vector<pop>::iterator p, q;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->popsize > threshold ){
			for(q=poplist.begin(); q != poplist.end(); q++){
				if(q->popsize > threshold && coexist(p->trait,q->trait) ){
					double sq = (1 - C(p->trait,q->trait)*K(p->trait)/K(q->trait) ), sp = (1 - C(p->trait,q->trait)*K(q->trait)/K(p->trait) ) ;
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


int finishline(vector<pop> &poplist, double line)
{
	int threshold = THRESHOLD;
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

void printtraits(vector<pop> &poplist)
{
	vector<pop>::iterator p;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->popsize > THRESHOLD) 
			printf("%lf ", p->trait);
	}
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


int traits_above_thresh(vector<pop> &poplist)
{
	vector<pop>::iterator p;
	int i = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->popsize > THRESHOLD)  i++;
	}
	return i;
}


double * gettraits(vector<pop> &poplist)
{
	double * dimorphism = (double *) 
		malloc((poplist.size())*sizeof(double));	
	vector<pop>::iterator p;
	int i = 0;
	for(p=poplist.begin(); p != poplist.end(); p++){
		if(p->popsize > THRESHOLD){
			 dimorphism[i] = p->trait;
			 i++;
		}
	}
	return dimorphism;
}



void boxplot(int x, int y, pngwriter &png)
{
	png.plot(x+1,y, 0., 0., 1.);
	png.plot(x-1,y, 0., 0., 1.);
	png.plot(x,y+1, 0., 0., 1.);
	png.plot(x,y-1, 0., 0., 1.);

	png.plot(x-1,y-1, 0., 0., 1.);
	png.plot(x+1,y+1, 0., 0., 1.);
	png.plot(x-1,y+1, 0., 0., 1.);
	png.plot(x+1,y-1, 0., 0., 1.);
/*
	png.plot(x+2,y, 0., 0., 1.);
	png.plot(x-2,y, 0., 0., 1.);
	png.plot(x,y+2, 0., 0., 1.);
	png.plot(x,y-2, 0., 0., 1.);
*/
}


void plotpng(vector<pop> &poplist, int time, pngwriter &png)
{
	vector<pop>::iterator p;
	int x=time, y;
	double ymax = .75;
	int shift = Y_PLOTSIZE/2;

	for(p=poplist.begin(); p != poplist.end(); p++){
		x = time;
		y = shift + (int) round((p->trait/ymax)*Y_PLOTSIZE/2);
		png.plot(x, y, 0.0, 0.0, 1.0);
		boxplot(x, y, png);
	}
	/* crude gridlines */
	png.plot(x, shift, .5, .5, .5);
	double secondline = .25;
	int grid = (int) ((secondline/ymax)*Y_PLOTSIZE/2);
	png.plot(x, shift+grid, .5, .5, .5);
	png.plot(x, shift-grid, .5, .5, .5);
}


