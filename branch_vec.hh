// Flags / Options
#define VERBOSE 1
#define CONTINUE 1
#define PNG_ON 0

// Parameters
#define MK (double) .5  // == 1/(2 sigma_k^2) ==> sigma_k^2 = 1
#define MC (double) 5. // == 1/(2 sigma_c^2) ==> sigma_c^2 = .2
#define SIGMA_MU (double) 0.05
#define Ko (int) 5000
#define MU (double) 1.e-4
#define Xo (double) -.5
#define SAMPLES (int) 1e4
#define MAXTIME 1e6
#define MAXTRIALS  1
#define THRESHOLD (int) Ko/10
#define X2 (double) -.5
#define N2o (double) 0.0
#define LINE (double) .6
#define EPSILON (double) 1e-9
#define R (double) 1.0
#define IK (double) 1.0/Ko
#define X_PLOTSIZE SAMPLES
#define Y_PLOTSIZE 500

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>

#include <iostream>
#include <list>
#include <vector>

#include <pngwriter.h>

using namespace std;

typedef struct {
		int    popsize;
		double trait;
		double birth;
		double death;
		double parent_trait;
		double compete;
}pop ;

typedef vector<double> CRow;

// Fast method functions
void print_C(vector<CRow> &cmatrix);
void grow_C(vector<pop> &poplist, vector<CRow> &cmatrix);
void shrink_C(vector<CRow> &cmatrix, int id);
void fast_update_rates(vector<pop> &poplist, vector<CRow> &cmatrix, int id, char type);

void event_and_rates(gsl_rng * rng, vector<pop> &poplist, double sum, vector<CRow> &cmatrix);
double init_compete(vector<pop> &poplist, double y);

// General methods
double C(double x, double y);
double K(double x);
double Kinv(double x);
double sumrates(vector<pop> &poplist);

double bdry(double x);
int coexist(double x, double y);
int branchcheck(vector<pop> &poplist, int threshold);
int finishline(vector<pop> &poplist, double line);
void averagelist(vector<pop> &poplist, double sampletime, double *mean, int samplenumber);

void printaverage(double *mean);
void printlist(vector<pop> &poplist, double sampletime);
void printfulllist(vector<pop> &poplist, double sampletime);
void printtraits(vector<pop> &poplist);
void printfreq(vector<pop> &poplist);

int mutual_invade(vector<pop> &poplist, int threshold);
int flagcheck(vector<pop> &poplist, double flag);
double flagget(vector<pop> &poplist, double flag);


// basicmethod functions
void update_rates(vector<pop> &poplist);
void event(gsl_rng * rng, vector<pop> &poplist, double sum);



void plotpng(vector<pop> &poplist, int time, pngwriter &png);

int traits_above_thresh(vector<pop> &poplist);


double * gettraits(vector<pop> &poplist);

int invade_pair(vector<pop> &poplist, int threshold, double * pair);
int branches(vector<pop> &poplist, int threshold, double * pair);

