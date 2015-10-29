#include "calc_xi_mm.h"

//This constant defines how much space the gsl workspace gets.
//It is repeated in another location in the code, so it is defined globally.
#define TOL 1e-6
#define workspace_size 8000
#define PI 3.141592653589793

//These are the parameters passed into the integrand.
//A spline and accelerator for interpolation and the radius we are evaluating xi at.
typedef struct integrand_params{
  gsl_spline *spline;
  gsl_interp_accel *acc;
  double r;
}integrand_params;

//These functions are a wrapper for the integral and the integrand itself.
double do_xi_mm_integral(double lkmin, double lkmax, integrand_params *params, gsl_integration_workspace *w);
double xi_mm_integrand(double lk, void * params);

//Here the xi_mm array is calculated.
int calc_xi_mm(double *radii,int numr,
	       double *z,int numz,double *k,int numk,
	       double *PK_in,int numPK1,int numPK2,
	       double **xi_mm){
  double r, xi_mm_at_r;//temporary variables
  int i,j,l;//iteration variables

  /*Note: PK comes in flattened, so we have to make a new variable
    for it in order to be able to index easily*/
  double PK[numPK1][numPK2];
  for(i=0;i<numPK1;i++)for(j=0;j<numPK2;j++)PK[i][j] = PK_in[i*numPK2 + j];

#pragma omp parallel shared(xi_mm,PK) private(i,j,l,r,xi_mm_at_r)
  {
  //Allocate memory for the spline and the accelerator
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_cspline,numk);
  gsl_interp_accel *acc
    = gsl_interp_accel_alloc();
  
  //Allocate workspace for the integration
  gsl_integration_workspace * workspace
    = gsl_integration_workspace_alloc(workspace_size);
  
  //Declare the function that we will integrate
  integrand_params*params = malloc(sizeof(integrand_params));
  params->acc=acc;

  //Find the log of the true kmin and kmax
  double lkmin,lkmax;
  lkmin = log(k[0]);
  lkmax = log(k[numk-1]);
  
  //Calculate xi_mm for eachredshift, mass, and radii
  for(i=0; i<numz; i++){
    //Create the spline for this redshift bin
    //Find the spline for the power spectrum. We are using a cubic spline
    gsl_spline_init(spline,k,PK[i],numk);
    params->spline=spline;
    
    //Loop over radii to find the correlation function there
#pragma omp for
    for(j=0;j<numr;j++){
      r = radii[j];
      params->r=r;
      
      xi_mm_at_r = do_xi_mm_integral(lkmin,lkmax,params,workspace);
      xi_mm[i][j] = xi_mm_at_r;
    }//end for j
  }//end for i

  //Free the gsl objects
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free (workspace);
  free(params);

  }//end parallel section

  return 0;
}

double do_xi_mm_integral(double lkmin, double lkmax, integrand_params *params,gsl_integration_workspace *w){
  //Define the function to be integrated.
  gsl_function F;
  F.function=&xi_mm_integrand;
  F.params=params;
  double result, abserr;
  //Integrate it.
  int status = gsl_integration_qag(&F,lkmin,lkmax,TOL,TOL/10.,workspace_size,6,w,&result,&abserr);
  //Check for errors.
  if(status){fprintf(stderr,"Error in gsl integration.\n");exit(status);}
  //printf("status = %d\nresult = %e\tabserr=%e\n",status,result,abserr);//Used for debugging.
  return result;
}

//double xi_mm_integrand(double lk, gsl_spline*spline, gsl_interp_accel *acc, double r){
double xi_mm_integrand(double lk, void *params){
  //Acquire the parameters
  integrand_params pars = *(integrand_params *)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double r = pars.r;

  //Return the evaluation
  double k = exp(lk);
  double x = k*r;
  double PK_k = gsl_spline_eval(spline,k,acc);
  return k*k*k*PK_k*sin(x)/x/(2.*PI*PI);
}
