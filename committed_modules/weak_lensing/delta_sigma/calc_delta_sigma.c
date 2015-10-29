#include "calc_delta_sigma.h"

//This constant defines how much space the gsl workspace gets.
//It is repeated in another location in the code, so it is defined globally.
#define workspace_size 8000
#define TOL 1e-3
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

//These are the parameters passed into the integrand.
//A spline and accelerator for interpolation and the radius we are evaluating xi at.
typedef struct integrand_params{
  gsl_spline *spline;
  gsl_interp_accel *acc;
  gsl_integration_workspace *w2;
  double rmax;
  double rmin;
  double r_p;
  double r_z;
  double R;
  double offset;
  double r_s,M,delta,H0,omega_m;//The cosmology and cluster features
}integrand_params;

//This is the analytic form of the sigma(R) 1halo value
double sigma_r_1halo_analytic(double r_s,double M,double omega_m,double H0,double delta,double R);

//These functions are wrappers for the integral and integrand istelf
double do_delta_sigma_integral(integrand_params*params,double rmin,
				    double rmax,double R,
				    gsl_integration_workspace*workspace1);
double integrand(double r, void * params);

//Here the surface density, \Delta\Sigma_z_M_R is calculated
int calc_delta_sigma(double M,double r_scale,double delta,
		     double omega_m,double H0,
		     double *z,int numz,
		     double *r,int numr,
		     double *sigma_r_in,double **delta_sigma){
  double R, s_lr,s_r;//Temp variables used in the integration
  int i,j,l;//iteration variables

  /*Note: sigma_r comes in flattened, so we have to make a new variable in order
    to be able to index it easily*/
  double sigma_r[numz][numr];
  for(i=0;i<numz;i++)
    for(j=0;j<numr;j++)
      sigma_r[i][j] = sigma_r_in[i*numr + j];
  
#pragma omp parallel shared(delta_sigma,sigma_r) private(i,j,l,R,s_lr,s_r)
  {
  //Allocate memory for the splines and the accelerators
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_cspline,numr);
  gsl_interp_accel *acc
    = gsl_interp_accel_alloc();
  
  //Allocate workspace for the integration
  //Note: we need two since we are doing a double integral
  gsl_integration_workspace * ws1
    = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace * ws2
    = gsl_integration_workspace_alloc(workspace_size);

  //Allocate space for the parameters for all integrands.
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline=spline;
  params->acc=acc;
  params->w2=ws2;

  //Put in the cosmology and the cluster features
  params->delta=delta;
  params->r_s=r_scale;
  params->M=M;
  params->omega_m=omega_m;
  params->H0=H0;
  
  //Find the min and max radii
  double rmin = r[0];
  double rmax = r[numr-1];
  params->rmax = rmax;
  params->rmin = rmin;
  
  //Calculate the surface density at each redshift, mass and radii
  for(i=0;i<numz;i++){
    //Create the spline for the current redshift and mass bin
    gsl_spline_init(spline,r,sigma_r[i],numr);
    
    //Loop over the radii
#pragma omp for 
    for(j=0;j<numr;j++){
      R = r[j];
      params->R=R;
      s_lr = do_delta_sigma_integral(params,rmin,rmax,R,ws1);
      s_r = sigma_r[i][j];
      delta_sigma[i][j] = s_lr-s_r;
      //The integral of the function has one less data point
      //Also, at very low R there can be an erroneous value
      if(delta_sigma[i][j]<0)delta_sigma[i][j] = 0;
    }//end for j
  }//end for i

  //Free the gsl objects
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free (ws1);
  gsl_integration_workspace_free (ws2);
  free(params);

  }/*end parallel section*/

  return 0;
}

double do_delta_sigma_integral(integrand_params*params,double rmin, 
			       double rmax,double R,
			       gsl_integration_workspace *workspace1){
  //Set the function we will integrate
  gsl_function F;
  F.function=&integrand;
  F.params=params;
  double result,abserr;
  //Integrate the function
  int status = gsl_integration_qag(&F,0,R,TOL,TOL/10.,workspace_size,
				   6,workspace1,&result,&abserr);
  if(status){fprintf(stderr,"Error after delta_sigma integration.\n");}
  result *= 2./(R*R);//the prefactors
  return result;
}

double integrand(double r_p,void*params){//The integral over r_p
  //Acquire the parameters
  integrand_params*pars = (integrand_params*)params;
  gsl_spline*spline = pars->spline;
  gsl_interp_accel*acc = pars->acc;
  double rmin = pars->rmin;
  double r_s = pars->r_s, M = pars->M, delta=pars->delta;
  double omega_m = pars->omega_m, H0 = pars->H0;
  double sigma_rp;
  if(r_p<rmin)//Check to see if we are outside of the interpolation
    sigma_rp = sigma_r_1halo_analytic(r_s,M,omega_m,H0,delta,r_p);
  else
    sigma_rp = gsl_spline_eval(spline,r_p,acc);
  //Return the evaluation
  return r_p*sigma_rp;
}

//The analytic form of sigma_r_1halo using NFW
double sigma_r_1halo_analytic(double r_s,double M,
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
    gx = (1 - 2./sqrt(1-x*x)*atanh(sqrt((1-x)/(1+x))))/(x*x-1);
  else if(x/r_s < 1e-6) // x is approximately equal to 1
    gx = 1./3.;
  else if(x>1)
    gx = (1 - 2./sqrt(x*x-1)* atan(sqrt((x-1)/(1+x))))/(x*x-1);
  //NOTE: an extra factor of 1e12 is divided in
  //this is because Sigma(r) has units of SM h/pc^2
  //So the Mpc^-2 must be converted to pc^-2
  return 2*r_s*delc*rhom*gx/(1e12);//SM h/pc^2
}
