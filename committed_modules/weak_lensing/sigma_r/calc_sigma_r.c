#include "calc_sigma_r.h"

//This constant defines how much space the gsl workspace gets.
//It is repeated in another location in the code, so it is defined globally.
#define TOL1 1e-3 //First level of the integral
#define TOL2 1e-6 //Second level of the integral for the offset
#define workspace_size 8000
#define PI 3.141592653589793

//These are physical constants
#define G 4.517e-48//Newton's G in Mpc^3/s^2/Solar Mass
#define Mpcperkm 3.241e-20//Mpc/km used to convert H0 to per seconds

//These are the parameters passed into the integrand.
//A spline and accelerator for interpolation and the radius we are evaluating xi at.
typedef struct integrand_params{
  gsl_spline *spline;
  gsl_interp_accel *acc;
  gsl_integration_workspace *w2;//WS for the inner integral
  double rmax;
  double rmin;
  double r_p;
  double r_z;
  double R;
  double offset;
  double r_s,M,delta,H0,omega_m;//The cosmology and cluster features
}integrand_params;

//This is the analytic form of the xi 1halo value
double xi_1halo_analytic(double r_s,double M,double omega_m,double H0,double delta,double R);

//These functions are wrappers for the integral and integrand istelf
double do_sigma_r_integral(integrand_params*params,double rmin,
			   double rmax,double R,double omega_m,double H0,
			   gsl_integration_workspace*workspace1);

double integrand(double rz, void * params);
double integrand_offset(double theta,void * params);

//Here the surface density, \Sigma_z_M_R is calculated
int calc_sigma_r(double M,double r_scale,double delta,
		 double *z,int numz,
		 double *r,int numr,double *xi_in,
		 double omega_m,double H0,double **sigma_r,
		 double offset){

  double R, s_r;//Temp variables used in the integration
  int i,j,l;//iteration variables
  
  /*Note: xi comes in flattened, so we have to make a new variable in order
    to be able to index it easily*/
  double xi[numz][numr];
  for(i=0;i<numz;i++)
    for(j=0;j<numr;j++)
      xi[i][j] = xi_in[i*numr + j];

#pragma omp parallel shared(sigma_r,xi) private(i,j,l,R,s_r)
  {
  //Allocate memory for the splines and the accelerators
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_cspline,numr);
  gsl_interp_accel *acc
    = gsl_interp_accel_alloc();
  
  //Allocate workspace for the integration
  /*Note: we need two since we are doing a double integral
    if there is an offset*/
  gsl_integration_workspace * ws1
    = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace * ws2
    = gsl_integration_workspace_alloc(workspace_size);

  //Allocate space for the parameters for all integrands.
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline=spline;
  params->acc=acc;
  params->w2=ws2;
  params->offset=offset;

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
    //Create the spline of xi for the current redshift
    gsl_spline_init(spline,r,xi[i],numr);
        
    //Loop over the radii
#pragma omp for 
    for(j=0;j<numr;j++){
      R = r[j];
      params->R=R;
      s_r = do_sigma_r_integral(params,rmin,rmax,R,omega_m,H0,ws1);
      if(s_r<0) s_r=0;
      sigma_r[i][j]=s_r;
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

double do_sigma_r_integral(integrand_params*params,double rmin,
			   double rmax,double R,double omega_m,double H0,
			   gsl_integration_workspace*workspace1){
  //Calculate the critical density in Solar Masses h^2/pc^2/Mpc
  //This is because the integration is over a domain with units of Mpc/h
  //however measurements of Sigma or Delta Sigma are in SM h/pc^2
  double rhocrit = 3.*(H0*H0*Mpcperkm*Mpcperkm)/(8.*PI*G)/(1e12)/(H0/100.*H0/100.);//SM h^2/pc^2/Mpc
  //Set the function we will integrate in part 2
  gsl_function F;
  F.function=&integrand;
  F.params=params;
  double result,abserr;
  //Integrate the function
  int status = gsl_integration_qag(&F,0,rmax/*sqrt(rmax*rmax-R*R)*/,TOL1,TOL1/10.,workspace_size,6,workspace1,&result,&abserr);
  if(status){fprintf(stderr,"Error after gsl integration.\n");}
  result *= rhocrit*omega_m*2.;//Symmetric about r_z=0, so factor of 2
  return result;
}

//This function is the integrand of Sigma (R)
double integrand(double rz, void * params){ 
  //Acquire the parameters
  integrand_params*pars = (integrand_params *)params;
  pars->r_z = rz;//Set the current r_z
  gsl_integration_workspace *w = pars->w2;
  gsl_function F;
  F.function = &integrand_offset;
  F.params=pars;
  double result,abserr;
  //int status = gsl_integration_qags(&F,0,2*PI,1e-6,1e-7,workspace_size,/*6,*/w,&result,&abserr);
  int status = gsl_integration_qag(&F,0,2*PI,TOL2,TOL2/10.,workspace_size,6,w,&result,&abserr);
  if(status){fprintf(stderr,"Error in offset integration.\n");exit(status);}
  return result/(2.*PI);
}

//This function is used when there is an offset
double integrand_offset(double theta, void * params){
  //Acquire the parameters
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;//Xi(R)
  gsl_interp_accel*acc = pars.acc;
  double R = pars.R;
  double rz = pars.r_z;
  double Rs = pars.offset;
  double rmax = pars.rmax;
  double rmin = pars.rmin;
  double r_s = pars.r_s, M = pars.M, delta=pars.delta;
  double omega_m = pars.omega_m, H0 = pars.H0;
  double arg = sqrt(rz*rz+R*R+Rs*Rs+2*R*Rs*cos(theta));
  //We have to account for going beyond the interpolation domain
  if(arg>rmax)
    return 0;//gsl_spline_eval(spline,rmax,acc);
  else if(arg<rmin)
    return xi_1halo_analytic(r_s,M,omega_m,H0,delta,arg);
  else
    return gsl_spline_eval(spline,arg,acc);
}

//The analytic form of xi_1halo using NFW
double xi_1halo_analytic(double r_s,double M,double omega_m,double H0,double delta,double R){
  double rhom = omega_m*3.*(H0*H0*Mpcperkm*Mpcperkm)/(8.*PI*G)
    /(H0/100.*H0/100.);//SM h^2/Mpc^3
  double rdelta = pow(M/(4./3.*PI*rhom*delta),1./3.);//Mpc/h
  double c = rdelta/r_s;//The concentration
  double fc = log(1+c)-c/(1+c);
  double x = R/r_s;
  return (M/(4.*PI*r_s*r_s*r_s*fc)/(x*(1+x)*(1+x)))/rhom - 1.0;
}
