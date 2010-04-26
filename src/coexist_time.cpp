#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "integrate.h"

	double Nhat(void * params){
		double k1 = ((double *) params)[0];
		double k2 = ((double *) params)[1];
		double c12  = ((double *) params)[2];
		//double Ko = ((double *) params)[3];
		return (k1+k2)/(1+c12);
	}
	double phat(void * params){
		double k1 = ((double *) params)[0];
		double k2 = ((double *) params)[1];
		double c12  = ((double *) params)[2];
		return (k1-c12*k2)/((1-c12)*(k1+k2));

	}

	double N (double x, void * params) {
		double k1 = ((double *) params)[0];
		double k2  = ((double *) params)[1];
		double c12  = ((double *) params)[2];

		double N = k1*k2 / 
				 (  k2*gsl_pow_2(x) 
				  + (k1+k2)*c12*x*(1-x) 
				  + k1*gsl_pow_2(1-x) ) ;
		return N;
	}
	double Np (double x, void * params) {
		double k1 = ((double *) params)[0];
		double k2  = ((double *) params)[1];
		double c12  = ((double *) params)[2];

		double Np = k1*k2 * 
				 pow((  k2*gsl_pow_2(x) 
				  + (k1+k2)*c12*x*(1-x) 
				  + k1*gsl_pow_2(1-x) ),-2) 
				  * ( k2*x 
				  + (k1+k2)*c12
				  - (k1+k2)*c12*x*2
				  - k1*2*(1-x) );
		return Np;
	}


	double A (double x, void * params) {
		double k1 = ((double *) params)[0];
		double c12  = ((double *) params)[2];

		double A = x*(1-N(x, params)*(x+c12*(1-x)) /k1) ;
//		double A = x*(1-Nhat(params)*(x+c12*(1-x)) /k1) ;
		return A;
	}

	double B (double x, void * params) {
		double k1 = ((double *) params)[0];
		double c12  = ((double *) params)[2];
		double Ko = ((double *) params)[3];

		double B = x*(1/N(x, params)+(x+c12*(1-x)) /k1)/Ko;
//		double B = x*(1/Nhat(params)+(x+c12*(1-x)) /k1)/Ko;
		return B;
	}
	double Bp (double x, void * params) {
		double k1 = ((double *) params)[0];
		double c12  = ((double *) params)[2];
		double Ko = ((double *) params)[3];
		double Bp = (1/(k1*Ko)) * (x*( pow(N(x,params),-2)*Np(x, params)+(1-c12))  +(1/N(x, params)+(x+c12*(1-x)) ) )  ;
//		double Bp = (1/(k1*Ko)) * (x*(1-c12)  +(1/Nhat(params)+(x+c12*(1-x)) ) )  ;
		return Bp;
	}

/*/
	double U (double x){return gsl_pow_4(x) - 4*gsl_pow_2(x) ; }
	double A (double x, void * params) { return -4*gsl_pow_3(x) + 8*x; }
	double B (double x, void * params) { return 2.0; }
	double UPsi (double x, void * params){ return exp( -U(x)/1.0 ); }
	double UPsi_inv (double x, void * params){ return exp( U(x)/1.0 ); }
*/


	double Psi_integrand (double x, void * params) { return 2*A(x, params)/B(x,params); }

	double Psi (double x, void * params) 
	{
		double a = ((double *) params)[4];
		return exp( integrate(&Psi_integrand, params, a, x) );
	}
	double Psi_inv (double x, void * params) 
	{
		double a = ((double *) params)[4];
		return exp(- integrate(&Psi_integrand, params, a, x) );
	}


	/* Psi(x)/B(x) */
	double Psi_B_integrand (double x, void * params) { return Psi(x, params)/B(x,params); }

	/* Int_a^x Psi(z)/B(z) dz */
	double Int_Psi_B (double x, void * params)
	{
		double a = ((double *) params)[4];
		return integrate(&Psi_B_integrand, params, a, x);
	}
	double Inner_integrand (double x, void * params) {	return Int_Psi_B(x,params)*Psi_inv(x,params); }


	double T(double x, void * params){
		double a = ((double *) params)[4];
		double b = ((double *) params)[5];
		return 2*(integrate(&Psi_inv, params,a,x)*integrate(&Inner_integrand, params,x,b) - integrate(&Psi_inv, params,x,b)*integrate(&Inner_integrand, params,a,x) )/integrate(&Psi_inv, params,a,b) ;
	}

	/* right boundary (b) aborbing, left (a) reflecting */
	double Tb(double x, void * params){
		double b = ((double *) params)[5];
		return 2*integrate(&Inner_integrand, params,x,b)  ;
	}


	/* Int_x^b Psi(z)/B(z) dz */
	double Int_Psi_B_b (double x, void * params)
	{
		double b = ((double *) params)[5];
		return integrate(&Psi_B_integrand, params, x, b);
	}
	double Inner_integrand_b (double x, void * params){	return Int_Psi_B_b(x,params)*Psi_inv(x,params); }


	/* left boundary (a) aborbing, right (b) reflecting */
	double Ta(double x, void * params){
		double a = ((double *) params)[4];
		return 2*integrate(&Inner_integrand_b, params,a,x)  ;
	}

	double prob_exit_a(double x, void * params){
		double a = ((double *) params)[4];
		double b = ((double *) params)[5];
		return integrate(&Psi, params, x, b)/integrate(&Psi, params, a, b);
	}


	double prob_exit_b(double x, void * params){
		double a = ((double *) params)[4];
		double b = ((double *) params)[5];
		return integrate(&Psi, params, a, x)/integrate(&Psi, params, a, b);
	}


	#define MK (double) .5  // == 1/(2 sigma_k^2) ==> sigma_k^2 = 1
	#define MC (double) 1. // == 1/(2 sigma_c^2) ==> sigma_c^2 = .2
	#define Ko (int) 500
	#define Xo (double) -.0059
	#define X2 -.0176

	double C(double x, double y){return exp( -gsl_pow_2(x-y)*MC );}
	double K(double x){return Ko*exp( -gsl_pow_2(x)*MK );}
	double S(double y, double x){ return 1 - K(x)*C(x,y)/K(y); }

	double energy(double x, void * params){
		return (A(x,params)-Bp(x,params)/4)/sqrt(B(x,params));
	}



	double bdry(double x){ 
		double sigma2c = 1/(2*MC), sigma2k = 1/(2*MK);
		return x*(1+sigma2c/sigma2k)/(1-sigma2c/sigma2k) ;
	}		
	double mirr(double x){ 
		double sigma2c = 1/(2*MC), sigma2k = 1/(2*MK);
		return x*(1-sigma2c/sigma2k)/(1+sigma2c/sigma2k) ;
	}

	int P2(double x1,double x2){ 
		return  ( (x1 > bdry(x2)) && (x1 > mirr(x2)) ) || ( (x1 < bdry(x2)) && (x1 < mirr(x2)) ) ; 
	}





void coexist_time (void)
{
	double k1 = K(Xo)/(double) Ko; 
	double k2 = K(X2)/(double) Ko; 
	double c12 = C(Xo, X2);

	double params[6] = {k1, k2, c12, (double) Ko, 1e-12, 1}, x = 1-1/(double) Ko;
//		double params2[6] = {K1, K2, C12, POPSIZE, 1e-6, 1}, x2 = 5/POPSIZE;
//		printf("%lf\n", k1/(c12*N(x,params) ) );

	double s =  (1 - (c12*N(x,params) )/k1 ) /
				(1 + (c12*N(x,params) )/k1 );


	double p_bar =	phat(params);
	double N_bar = Ko*Nhat(params); 
	double N_1 = Ko*(k1-c12*k2)/(1-c12*c12);
	double N_2 = Ko*(k2-c12*k1)/(1-c12*c12);

//		printf("#  p_bar = %lf, N_bar = %lf, N_1 = %lf, N_2 = %lf\n", p_bar, N_bar, N_1, N_2);
//		printf("#  s = %lf S12 = %lf, S21 = %lf\n", s, S(Xo,X2), S(X2,Xo) );
	printf("#  X1 = %lf X2 = %lf, MC = %lf, MK = %lf, Ko = %d\n", Xo, X2, MC, MK, Ko);
//		printf("#  k1 = %lf k2 = %lf, c12 = %lf, x =%lf\n", k1, k2, c12, x);

//		printf("#  T pop1 = %lf, %lf, %lf, %lf\n",  T(.5, params), params[0], params[1], params[2] );

/*

//		printf("%lf, %lf\n", exp( energy(1,params)-energy(p_bar,params) ), energy(1,params));
	printf("#  pi_a = %lf\n",  prob_exit_a(x, params) );
	printf("#  pi_b = %lf\n",  prob_exit_b(x, params) );
	printf("#  T pop1 = %lf\n",  T(x, params) );
//		printf("#  T pop2 = %lf\n",  T(x2, params2) );
	printf("#  Ta = %lf\n",  Ta(x, params) );
	printf("#  Tb = %lf\n",  Tb(x, params) );
	printf("#  Ta from pbar = %lf\n",  Ta(p_bar, params) );
	printf("#  Tb from pbar = %lf\n",  Tb(p_bar, params) );
	printf("#  1-exp(-2s)*Ta(pbar) = %lf\n",  (1-exp(-2*s) ) * Ta(p_bar, params) );
	printf("#  kimura u = 1-exp(-2s) = %lf\n",  (1-exp(-2*s) ) );
	printf("#  1/(1/Ta + 1/Tb)  = %lf\n",   1/(1/Ta(x, params)+1/Tb(x, params) ));
	printf("#  1/(1/Ta(pbar) + 1/Tb(pbar))  = %lf\n",   1/(1/Ta(p_bar, params)+1/Tb(p_bar, params) ));
	printf("#  u* harm. mean from pbar = %lf\n",  (1-exp(-2*s) ) * 1/(1/Ta(p_bar, params)+1/Tb(p_bar, params) ));
	printf("#  u*harm. mean from x = %lf\n",  (1-exp(-2*s) ) * 1/(1/Ta(x, params)+1/Tb(x, params) ));
	printf("#  weighted mean = %lf\n",  1/(1/((1-exp(-2*s))*Ta(x, params))+1/(exp(-2*s)*Tb(x, params))) );

	printf("#  pi weights = %lf\n",  1/(1/(Ta(x,params)*prob_exit_a(x, params)) + 1/(Tb(x,params)*prob_exit_b(x, params)) ) );

	double p = 0.01;
	printf("#  A = %lf, B = %lf, Bp = %lf\n", A(p, params), B(p, params), Bp(p, params));

*/		
	//printf("Arrhenius = %lf\n",  2*M_PI*exp( (U(0) -U(-sqrt(2)) )/1  )*(1/sqrt(16))*(1/sqrt(8)) );
	//printf("SD pop1 = %lf\n",  sqrt( SecondMoment(x,params) - gsl_pow_2( T(x, params) ) ));
	
	

	int i,j;
	double x1=.06, x2= .02;

	params[0] = K(x1)/ (double) Ko; 
	params[1] = K(x2)/ (double) Ko; 
	params[2] = C(x1,x2); 

	printf("phat %lf\n", phat(params));
//		printf("%d %d\n", P2(.06, .02), P2(-.1, -.1) );
//		for(i=0;i<20;i++) printf("%lf\n", bdry(.02*i+-.2) );

	x1=-.2;
	double ph;
	for(i=0; i<20; i++){
		x1 += .02;
		x2 = -.2;
		params[0] = K(x1)/ (double) Ko; 
		for(j=0; j<20; j++){
			x2 += .02;
			if( P2(x1,x2) ){
				params[1] = K(x2)/ (double) Ko; 
				params[2] = C(x1,x2);
				ph = GSL_MAX(phat(params), 1/(double) Ko);
				printf("%lf, %lf, %lf %lf\n", T(ph, params), x1, x2, ph);
			}
			else printf("NA, %lf, %lf\n", x1, x2 );
		}
	}
}	

	/* 
	 *		U(x) = x^4 - 4 x^2
	 *		U'(x) = 4x^3 - 8 x
	 *		min's at +/- sqrt(2), max at 0
	 *		U(0) = 0, U(sqrt(2)) = -4
	 *
	 *		\partial_x [U'(x) p(x,t)] + D \partial_x^2 p(x,t)
	 *
	 *		Psi(x) = exp(-U(x)/D )
	 *
	 *		Start at -sqrt(2), get to max: 0
	 *		T ~ 2 alpha beta pi exp( ( U(0) - U(- sqrt(2)) )/D )
	 *		1/alpha^2 = U''(sqrt(2) ) = 12*sqrt(2)^2-8 = 16
	 *		1/beta^2 =  U''( 0 ) = -8 
	 *
	 *
	 * */

/* ---------------------------- */
/* Psi(x) T(x) /B(x) */
double Psi_T_B_integrand (double x, void * params) { return T(x,params)*Psi(x, params)/B(x,params); }

/* Int_a^x Psi(z)/B(z) dz */
double Int_Psi_T_B (double x, void * params)
{
	double a = ((double *) params)[4];
	return integrate(&Psi_T_B_integrand, params, a, x);
}
double SecondMoment_integrand (double x, void * params) {	return Int_Psi_T_B(x,params)*Psi_inv(x,params); }

double SecondMoment(double x, void * params){
	double a = ((double *) params)[4];
	double b = ((double *) params)[5];
	return 2*(integrate(&Psi_inv, params,a,x)*integrate(&SecondMoment_integrand, params,x,b) 
	- integrate(&Psi_inv, params,x,b)*integrate(&SecondMoment_integrand, params,a,x) )/integrate(&Psi_inv, params,a,b) ;
}



//Vp
//-(((4*c-4)*k+(4*c-4)*j)*r*gsl_pow(x,3)+((11-7*c)*k+(-3*c-1)*j)*r*gsl_pow(x,2)+((2*c-8)*k+(2-4*c)*j)*r*x+(k+c*j)*r)/(4*j*k*sqrt(-(((c-1)*k+(c-1)*j)*r*gsl_pow(x,3)+((2-c)*k-j)*r*gsl_pow(x,2)+(-k-c*j)*r*x)/(j*k*K))*K)

