#include "calc_ave_delta_sigma.h"

//This constant defines how much space the gsl workspace gets.
//It is repeated in another location in the code, so it is defined globally.
#define workspace_size 8000
#define TOL 1e-2
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

//These are the parameters passed into the integrand.
//A spline and accelerator for interpolation and the radius we are evaluating delta_sigma at.
typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace*w;
  double bin_max;
  double bin_min;
}integrand_params;

//These functions are wrappers for the integral and integrand istelf
double do_ave_delsig_integral(integrand_params*params);
double integrand(double r, void * params);

//Here the average delta sigma (z,M,R) is calculated
int calc_ave_delsig(double *z,int numz,
		    double *r,int numr,
		    double *bin_edges,int num_bins,
		    double *delta_sigma_in,
		    double **ave_delsig){
  double bin_min,bin_max, ads; //Temporary variables
  int i,j,l; //iteration variables
  
  /*Note: delta_sigma comes in flattened, 
    so we have to make a new variable in order 
    to be able to index it easily*/
  double delta_sigma[numz][numr];
  for(i=0;i<numz;i++)
    for(j=0;j<numr;j++)
      delta_sigma[i][j] = delta_sigma_in[i*numr + j];

#pragma omp parallel shared(ave_delsig,delta_sigma) private(i,j,l,bin_min,bin_max,ads)
  {
  //Allocate memory for the splines and the accelerators
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_cspline,numr);
  gsl_interp_accel *acc
    = gsl_interp_accel_alloc();

  //Allocate workspace for the integration
  gsl_integration_workspace * ws
    = gsl_integration_workspace_alloc(workspace_size);

  //Allocate space for the parameters for all integrands.
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline=spline;
  params->acc=acc;
  params->w=ws;

  //Calculate the average surface density at each redshift, mass and radii
  for(i=0;i<numz;i++){
    //Create the spline for the current redshift and mass bin
    gsl_spline_init(spline,r,delta_sigma[i],numr);

    //Loop over the radii
#pragma omp for
    for(j=0;j<num_bins;j++){
      bin_min = r[j];
      bin_max = r[j+1];
      params->bin_min=bin_min;
      params->bin_max=bin_max;
      ads = do_ave_delsig_integral(params);
      if(ads<0)ads=0;
      ave_delsig[i][j]=ads;
    }//end for j
  }//end for i

  //Free the gsl objects
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free (ws);
  free(params);

  }//End parallel section
  return 0;
}

//Perform the integral over Delta Sigma to calculate the 
//average DeltaSigma in a radial bin
double do_ave_delsig_integral(integrand_params*params){
  gsl_function F;
  F.function = &integrand;
  F.params = params;
  double bin_min = params->bin_min;
  double bin_max = params->bin_max;
  double result=0,abserr=0;
  int status = gsl_integration_qag(&F,bin_min,bin_max,TOL,TOL/10.,workspace_size,
				   6,params->w,&result,&abserr);
  if(status){fprintf(stderr,"Error after ave_delsig integration.\n");}
  result *=2./(bin_max*bin_max-bin_min*bin_min);
  return result;
}

//The integrand is just R*DeltaSigma(R)
double integrand(double r,void*params){
  integrand_params*pars=(integrand_params*)params;
  gsl_spline*spline = pars->spline;
  gsl_interp_accel*acc = pars->acc;
  double delta_sigma;
  delta_sigma = gsl_spline_eval(spline,r,acc);
  return r*delta_sigma;
}
/*
//This is the analytic form of the Delta Sigma(R) 1halo value
double delta_sigma_r_1halo_analytic(double r_s,double M,
double omega_m,double H0,
double delta,double R){
double rhom = omega_m*3.*(H0*H0*Mpcperkm*Mpcperkm)/(8.*PI*G)
/(H0/100.*H0/100.);//SM h^2/Mpc^3
double rdelta = pow(M/(4./3.*PI*rhom*delta),1./3.);//Mpc/h
double c = rdelta/r_s;//The concentration
double delc = (delta/3.)*c*c*c/(log(1+c)-c/(1+c));
double x = R/r_s;
double gx = 0;
if(x<1)
gx = 8.*atanh(sqrt((1-x)/(1+x)))/(x*x*sqrt(1-x*x)) 
+ 4.*log(x/2.)/(x*x) - 2./(x*x-1)
+ 4.*atanh(sqrt((1-x)/(1+x)))/((x*x-1)*sqrt(1-x*x));
else if(x/r_s < 1e-6) // x is approximately equal to 1
gx = 10./3. + 4.*log(1./2.);
else if(x>1)
gx = 8.*atan(sqrt((x-1)/(1+x)))/(x*x*sqrt(x*x-1)) 
+ 4.*log(x/2.)/(x*x) - 2./(x*x-1)
+ 4.*atan(sqrt((x-1)/(1+x)))/((x*x-1)*sqrt(x*x-1));
//NOTE: an extra factor of 1e12 is divided in
//this is because Delta Sigma(r) has units of SM h/pc^2
//So the Mpc^-2 must be converted to pc^-2
return r_s*delc*rhom*gx/(1e12);//SM h/pc^2
}
*/
