#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "branch_vec.hh"


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


#define GRID 20	

	int analytic_contours (par_list * pars, double *xval, double *yval, double *times)
	{
		double ko = pars->ko, x1, x2, ph;
		int i, j;
		/* hitting time parameters 	k1,k2,c12,Ko,left_bdry,right_bdry*/
		double params[6] = {0, 0, 0, ko, 1e-12, 1};

		x1= - pars->xo;
		#pragma omp for private(x1, x2)
		for(i=0; i<GRID; i++){
			x1 += 2*pars->xo/GRID;
			x2 = -pars->xo;
			params[0] = K(x1, pars)/ ko; 
			for(j=0; j<GRID; j++){
				x2 += 2*pars->xo/GRID;
				if( coexist(x1,x2, pars) ){
					params[1] = K(x2, pars)/ ko; 
					params[2] = C(x1,x2, pars);
					ph = GSL_MAX(phat(params), 1/ko);
					double expected_coexistence = T(ph, params);

					xval[GRID*i+j] = x1;
					yval[GRID*i+j] = x2;
					times[GRID*i+j] = expected_coexistence;

//					printf("%lf, %lf, %lf %e\n", expected_coexistence, x1, x2, ph);
				}
				else{
					xval[GRID*i+j] = x1;
					yval[GRID*i+j] = x2;
					times[GRID*i+j] = 0;
//					printf("0, %lf, %lf\n", x1, x2);
				}
			}
		}
		return 0;
	}	


void analytic_contours_wrapper(double *sigma_mu, double *mu, double *sigma_c2, double *sigma_k2, double *ko, double *xo, double * times, double * xvals, double * yvals)
{
	double mc = 1 / (2 * *sigma_c2);
	double mk = 1 / (2 * *sigma_k2);
	par_list p = {*sigma_mu, *mu, mc, mk, *ko, 1 / *ko, *xo, NULL};
	par_list * pars = &p;

	analytic_contours(pars, xvals, yvals, times);

}

