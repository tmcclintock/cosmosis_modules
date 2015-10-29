/*
  This module will calculate the tinker bias.
 */
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_tinker_bias.h"

//These are some section names we will use
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * dist = DISTANCES_SECTION;
const char * setup_name = "setup_delta_sigma";
const char * matter_power_blockname = "matter_power_lin";

/* This struct will hold the data from the options.*/
typedef struct tinker_bias_config{
  double delta;//The overdensity value
} tinker_bias_config;

/* In setup, the options are obtained and passed 
   to execute via the config object.*/
void * setup(c_datablock * options){
  DATABLOCK_STATUS status=0;
  tinker_bias_config * config = malloc(sizeof(tinker_bias_config));
  status |= c_datablock_get_double(options,OPTION_SECTION,
				   "delta",&(config->delta));
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
  DATABLOCK_STATUS status=0;//used to track the status of the program
  int extents[2];//used to read in 2D arrays
  int ndims;//used to read in 2D arrays

  //Obtain the options from config.
  tinker_bias_config * config = (tinker_bias_config *)config_in;
  double delta = config->delta;

  //Acquire the mass
  double M;
  status |= c_datablock_get_double(block,cosmo,"log10_cluster_mass",&M);
  M = pow(10,M);

  //Acquire the redshift samples and the matter power spectrum
  double *z;
  double *k;
  int numz, numk;
  status |= c_datablock_get_double_array_1d(block,matter_power_blockname,"z",&z,&numz);
  status |= c_datablock_get_double_array_1d(block,matter_power_blockname,"k_h",&k,&numk);
  status |= c_datablock_get_array_ndim(block,matter_power_blockname,"p_k",&ndims);
  status |= c_datablock_get_double_array_shape(block,matter_power_blockname,"p_k",ndims,extents);
  if(status){fprintf(stderr,"Error on reading in PK dimensions.\n");exit(status);}
  int numPK1 = extents[0],numPK2=extents[1];
  double PK[numPK1][numPK2];
  status |= c_datablock_get_double_array(block,matter_power_blockname,
					 "p_k",(double *)PK,ndims,extents);
  if(status){fprintf(stderr,"Error on reading in PK.\n");exit(status);}
  if(numz!=numPK1||numk!=numPK2){
    fprintf(stderr,"Error, numz!=numPk1 or numk!=numPK2.\n");
    fprintf(stderr,"numz=%d\tnumPK1=%d\tnumk=%d\tnumPK2=%d\n",numz,numPK1,numk,numPK2);
    exit(1);
  }//Check sizes of arrays
  
  //Acquire the cosmological parameters
  double omega_m, H0;
  status |= c_datablock_get_double(block,cosmo,"omega_m",&omega_m);
  status |= c_datablock_get_double(block,cosmo,"hubble",&H0);
  if(status){fprintf(stderr,"Error on reading in cosmology.\n");exit(status);}

  //Declare and allocate the final array to be written to the block
  int ndimbias = 1;
  int extentsbias[] = {numz};
  double *bias = (double *)malloc(numz*sizeof(double));
  double *sigma2m= (double *)malloc(numz*sizeof(double));
  
  //Do the bias calculation
  //Note: we have to flatten the PK array
  status |= calc_tinker_bias(delta,M,z,numz,k,numk,(double*)PK,H0,omega_m,bias,sigma2m);

  //Write the bias to the block
  char *section="tinker_bias";char *name="bias"; char *signame="sigma2m";
  status |= c_datablock_put_double_array(block,section,name,
					 bias,ndimbias,extentsbias);
  status |= c_datablock_put_double_array(block,section,//Has the same dimenstions as bias
					 signame,sigma2m,ndimbias,extentsbias);

  //Free the bias array
  free(bias);free(sigma2m);

  return 0;//No likelihood from this module
}

//In cleanup we have to free the configuration data.
int cleanup(void * config_in){
  tinker_bias_config * config = (tinker_bias_config *)config_in;
  free(config);
  return 0;
}
