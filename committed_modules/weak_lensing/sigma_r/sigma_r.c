/*
  This module will calculate the surface density, delta sigma
 */
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_sigma_r.h"

//These are some section names we will use
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * dist = DISTANCES_SECTION;
const char * setup_name = "setup_delta_sigma";//Hard coded, for now
const char * xi_hm = "xi_hm";

/* This struct will hold the data from the options.*/
typedef struct sigma_r_config{
  double offset;
  double delta;
} sigma_r_config;

/* In setup, the options are obtained and passed to execute via the config object.*/
void * setup(c_datablock * options){
  DATABLOCK_STATUS status=0;
  sigma_r_config * config = malloc(sizeof(sigma_r_config));
  status |= c_datablock_get_double(options,OPTION_SECTION,"offset",&(config->offset));
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
  sigma_r_config * config = (sigma_r_config *)config_in;
  double offset = config->offset;
  double delta = config->delta;

  //Acquire the cluster mass and scale radius
  double M;
  double r_scale;
  status |= c_datablock_get_double(block,cosmo,"log10_cluster_mass",&M);
  status |= c_datablock_get_double(block,cosmo,"scale_radius",&r_scale);
  M = pow(10,M);

  //Acquire the radii and redshift to evaluate the correlation function at.
  double *r,*z;
  int numr;
  int numz;
  status |= c_datablock_get_double_array_1d(block,dist,"z",&z,&numz);
  status |= c_datablock_get_double_array_1d(block,setup_name,"radii",&r,&numr);

  //Acquire omega_m, the only cosmological parameter that comes into play here
  double omega_m, H0;
  status |= c_datablock_get_double(block,cosmo,"omega_m",&omega_m);
  status |= c_datablock_get_double(block,cosmo,"hubble",&H0);
  
  //Acquire the correlation function, xi, from the block
  status |= c_datablock_get_array_ndim(block,xi_hm,"xi_hm",&ndims);
  status |= c_datablock_get_double_array_shape(block,xi_hm,"xi_hm",ndims,extents);
  if(status){fprintf(stderr,"Error on reading in xi (2D) dimensions.\n");exit(status);}
  int numxi0 = extents[0],numxi1=extents[1];
  double xi[numxi0][numxi1];
  status |= c_datablock_get_double_array(block,xi_hm,"xi_hm",(double *)xi,ndims,extents);
  if(status){fprintf(stderr,"Error on reading in xi.\n");exit(status);}


  //Declare and allocate the resulting array for the surface density, \bar{Delta Sigma}(<R),
  //and \bar{Delta Sigma}(R)
  int ndimsd=2;
  int extentssd[]={numz,numr};
  double **sigma_r;

  sigma_r = (double **)malloc(numz*sizeof(double*));
  for(i=0;i<numz;i++)
    sigma_r[i] = (double *)malloc(numr*sizeof(double));
  
  //Do the surface density calculation
  //Note: we have to flatten the xi array
  status |= calc_sigma_r(M,r_scale,delta,z,numz,r,numr,(double *)xi,omega_m,H0,sigma_r,offset);

  //We have to flatten the sigma_r array in order to write to the block
  double sigma_r_out[numz][numr];
  for(i=0;i<numz;i++)
    for(j=0;j<numr;j++)
      sigma_r_out[i][j]=sigma_r[i][j];
  
  //Write sigma_r to the block
  char *name = "sigma_r";
  status |= c_datablock_put_double_array(block,name,name,(double *)sigma_r_out,ndimsd,extentssd);
  if(status){fprintf(stderr,"Error on writing out sigma_r to the block.\n");exit(status);}

  //Free the sigma_r array
  for(i=0;i<numz;i++)
    free(sigma_r[i]);
  free(sigma_r);

  return 0;//No likelihood from this module
}

//In cleanup we have to free the configuration data.
int cleanup(void * config_in){
  sigma_r_config * config = (sigma_r_config *)config_in;
  free(config);
  return 0;
}
