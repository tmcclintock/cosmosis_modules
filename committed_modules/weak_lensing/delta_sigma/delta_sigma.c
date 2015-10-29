/*
  This module will calculate the surface density, delta sigma
 */
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_delta_sigma.h"

//These are some section names we will use
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * dist = DISTANCES_SECTION;
const char * setup_name = "setup_delta_sigma";//Hard coded, for now
const char * sigma_r_name = "sigma_r";

/* This struct will hold the data from the options.*/
typedef struct delta_sigma_config{
  double delta;
} delta_sigma_config;

/* In setup, the options are obtained and passed to execute via the config object.*/
void * setup(c_datablock * options){
  DATABLOCK_STATUS status=0;
  delta_sigma_config*config = malloc(sizeof(delta_sigma_config));
  status |= c_datablock_get_double(options,OPTION_SECTION,"delta",&(config->delta));
  if(status){fprintf(stderr,"Error in setup.\n");exit(status);}
  return config;
}

/* In execute, the options are taken in.
   Next, data is read from the block that was
   generated previously.
   Then, the calculation is kicked off.
   Then, the program writes to the block.*/
int execute(c_datablock * block, void * config_in){
  int i,j,l;//iteration variables
  DATABLOCK_STATUS status=0;
  int extents[]={0,0,0};//Used to read in >1D arrays. Note: the extents are zerod out in order to test array lengths
  int ndims;//Used to read in >1D arrays

  //Obtain the options from config.
  delta_sigma_config*config = (delta_sigma_config*)config_in;
  double delta = config->delta;//The overdensity factor (usually 200)

  //Acquire the cluster mass and scale radius
  double M, r_scale;
  status |= c_datablock_get_double(block,cosmo,"log10_cluster_mass",&M);
  status |= c_datablock_get_double(block,cosmo,"scale_radius",&r_scale);
  M = pow(10,M);

  //Acquire the radii to evaluate sigma_r
  double *r;
  int numr;
  status |= c_datablock_get_double_array_1d(block,setup_name,"radii",&r,&numr);

  //Acquire omega_m and H0
  double omega_m, H0;
  status |= c_datablock_get_double(block,cosmo,"omega_m",&omega_m);
  status |= c_datablock_get_double(block,cosmo,"hubble",&H0);

  //Acquire the redshift samples and the matter density
  double *z;
  int numz;
  status |= c_datablock_get_double_array_1d(block,dist,"z",&z,&numz);
  
  //Acquire the surface density, sigma_r, from the block
  status |= c_datablock_get_array_ndim(block,sigma_r_name,"sigma_r",&ndims);
  status |= c_datablock_get_double_array_shape(block,sigma_r_name,"sigma_r",ndims,extents);
  if(status){fprintf(stderr,"Error on reading in sigma_r (2D) dimensions.\n");exit(status);}
  int numsigr0 = extents[0],numsigr1=extents[1];
  double sigma_r[numsigr0][numsigr1];
  status |= c_datablock_get_double_array(block,sigma_r_name,"sigma_r",(double *)sigma_r,ndims,extents);
  if(status){fprintf(stderr,"Error on reading in sigma_r.\n");exit(status);}

  //Declare and allocate the resulting array for the surface density
  int ndimsd=2;
  int extentssd[]={numz,numr};
  double **delta_sigma;
  delta_sigma = (double **)malloc(numz*sizeof(double*));
  for(i=0;i<numz;i++)
    delta_sigma[i] = (double *)malloc(numr*sizeof(double));
  
  //Do the surface overdensity calculation
  //Note: we have to flatten the sigma_r array
  status |= calc_delta_sigma(M,r_scale,delta,omega_m,H0,z,numz,r,numr,(double *)sigma_r,delta_sigma);

  //We have to flatten the delta_sigma array in order to write to the block
  double delta_sigma_out[numz][numr];
  for(i=0;i<numz;i++)
    for(j=0;j<numr;j++)
      delta_sigma_out[i][j]=delta_sigma[i][j];
  
  //Write delta_sigma to the block
  char *name = "delta_sigma";
  status |= c_datablock_put_double_array(block,name,name,(double *)delta_sigma_out,ndimsd,extentssd);
  if(status){fprintf(stderr,"Error on writing out delta_sigma to the block.\n");exit(status);}

  //Free the delta_sigma array
  for(i=0;i<numz;i++)
    free(delta_sigma[i]);
  free(delta_sigma);

  return 0;//No likelihood from this module
}

//In cleanup we have to free the configuration data.
int cleanup(void * config_in){
  delta_sigma_config * config = (delta_sigma_config *)config_in;
  free(config);
  return 0;
}
