#include "calc_tinker_bias.h"

//This constant defines how much space the gsl workspace gets.
//It is repeated in another location in the code, so it is defined globally.
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
  double R;
}integrand_params;

//These are wrappers for internal functions
int calc_nu(double M,double*z,int numz,double*k,int numk,
	    double*PK_in,double H0,double omega_m,double*nu,
	    double*sigma2m);
double do_sigma2_integral(double lkmin,double lkmax,
			  integrand_params*params,
			  gsl_integration_workspace*workspace);
double sigma2_integrand(double lk,void*params);

//Here the tinker bias is calculated
int calc_tinker_bias(double delta,double M,double*z,int numz,
		     double*k,int numk,double*PK_in,
		     double H0,double omega_m,double*bias,
		     double*sigma2m){  
  int i,j,l;//iteration variables
  int status=0;//an error checking variable

  //Calculate the variables associated with delta
  double delta_c = 1.686;//Critical collapse density
  double y = log10(delta);
  double xp = exp(-1.0*pow(4./y,4.));
  double A = 1.+0.24*y*xp, a = 0.44*y-0.88;
  double B = 0.183, b = 1.5;
  double C = 0.019+0.107*y+0.19*xp, c = 2.4;

  //Calculate the linear density field peak
  double *nu=(double*)malloc(numz*sizeof(double));
  status |= calc_nu(M,z,numz,k,numk,PK_in,H0,omega_m,nu,sigma2m);

  //Calculate the bias function
  for(i=0;i<numz;i++){
    bias[i] = 1 
      - A*pow(nu[i],a) / (pow(nu[i],a)+pow(delta_c,a))
      + B*pow(nu[i],b) 
      + C*pow(nu[i],c);
  }//end i
  
  //Free nu
  free(nu);

  return 0;
}

//This function calculates the peak of the linear density field
int calc_nu(double M,double*z,int numz,double*k,int numk,
	    double*PK_in,double H0,double omega_m,double*nu,double*sigma2m){
  int i,j,l;//iteration variables
  int status=0;//an error checking variable
  
  /*Note: PK came in flattened, so we have to make a new variable
    for it in order to be able to index it easily*/
  double PK[numz][numk];
  for(i=0;i<numz;i++)for(j=0;j<numk;j++)PK[i][j] = PK_in[i*numk + j];

  double delta_c = 1.686;//Critical collapse density

  //Calculate the lagrangian radius for each mass
  double R;
  double rhocrit=3.*(H0*Mpcperkm*H0*Mpcperkm)/(8.*PI*G)
    /(H0/100.*H0/100.);//SM h^2/Mpc^3
  double rhom = omega_m*rhocrit;//Comoving rho_matter
  R=pow(M/(4./3.*PI*rhom),1./3.);//Mpc/h

  //Calculate the minimum and maximum k values
  double lkmin = log(k[0]);
  double lkmax = log(k[numk-1]);

  //Declare the variance of the linear matter power
  double sigma2=0;//is a temporary variable

  //Allocate memory for the spline and the accelerator
  gsl_spline *spline
    = gsl_spline_alloc(gsl_interp_cspline,numk);
  gsl_interp_accel *acc
    = gsl_interp_accel_alloc();

  //Allocate workspace for the integration
  gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(workspace_size);

  //Declare the function that we will integrate
  integrand_params *params = (integrand_params*)malloc(sizeof(integrand_params));
  params->spline=spline;
  params->acc=acc;
  params->R=R;

  //For each redshift calculate the peak
  for(i=0;i<numz;i++){
    //Create the spline for this redshift bin
    //Find the spline for the power spectrum. We are using a cubic spline
    gsl_spline_init(spline,k,PK[i],numk);
    
    sigma2 = do_sigma2_integral(lkmin,lkmax,params,workspace);
    sigma2m[i] = sigma2;
    nu[i] = delta_c/sqrt(sigma2);
  }//end i

  //Free the gsl objects
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free (workspace);
  free(params);

  return 0;
}

//Here, the integral for sigma2 is actually performed
double do_sigma2_integral(double lkmin,double lkmax,
			  integrand_params*params,
			  gsl_integration_workspace*workspace){
  //Define the function to be integrated
  gsl_function F;
  F.function = &sigma2_integrand;
  F.params = params;
  double result, abserr;
  //Integrate it
  int status = gsl_integration_qag(&F,lkmin,lkmax,1e-9,1e-10,
				    workspace_size,6,workspace,
				   &result,&abserr);
  if(status){fprintf(stderr,"Error in gsl integration.\n");exit(status);}
  return result;
}

double sigma2_integrand(double lk,void*params){
  //Acquire the parameters
  integrand_params pars = *(integrand_params *)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double R = pars.R;

  double k = exp(lk);
  double x = k*R;
  double PK_k = gsl_spline_eval(spline,k,acc);
  double w = 3.0/x/x/x*(sin(x)-x*cos(x));//Window function
  return k*k*k*w*w*PK_k/(2.*PI*PI);
}
