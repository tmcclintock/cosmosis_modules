/*
  This module will calculate the average surface density, delta sigma (R)
 */
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_ave_delta_sigma.h"

//These are some section names we will use
const char * dist = DISTANCES_SECTION;
const char * setup_name = "setup_delta_sigma";//Hard coded, for now
const char * delta_sigma_name = "delta_sigma";

/* This struct will hold the data from the options.*/
typedef struct ave_delsig_config{
  double bin_min, bin_max;
  int num_bins;
} ave_delsig_config;

/* In setup, the options are obtained and passed to execute via the config object.*/
void * setup(c_datablock * options){
  DATABLOCK_STATUS status=0;
  ave_delsig_config*config = malloc(sizeof(ave_delsig_config));
  status |= c_datablock_get_double(options,OPTION_SECTION,"bin_min",&(config->bin_min));
  status |= c_datablock_get_double(options,OPTION_SECTION,"bin_max",&(config->bin_max));
  status |= c_datablock_get_int(options,OPTION_SECTION,"num_bins",&(config->num_bins));
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
  int extents[]={0,0};//Used to read in >1D arrays. Note: the extents are zerod out in order to test array lengths
  int ndims;//Used to read in >1D arrays

  //Obtain the options from config.
  ave_delsig_config*config = (ave_delsig_config*)config_in;
  double bin_min = config->bin_min;
  double bin_max = config->bin_max;
  int num_bins = config->num_bins;

  //Calculate the bin edges
  double bin_edges[num_bins+1];
  double dlr = (log(bin_max)-log(bin_min))/num_bins;
  for(i=0;i<num_bins+1;i++)
    bin_edges[i] = exp(log(bin_min)+i*dlr);

  //Acquire the redshift samples and the matter density
  double *z;
  int numz;
  status |= c_datablock_get_double_array_1d(block,dist,"z",&z,&numz);

  //Acquire the radii to evaluate delta_sigma
  double *r;
  int numr;
  status |= c_datablock_get_double_array_1d(block,setup_name,"radii",&r,&numr);

  //Check that the minimum and maximum bins are within the DeltaSigma radii
  if(bin_min < r[0]){fprintf(stderr,"Error: bin_min < r[0].\n");exit(2);}
  if(bin_max < r[num_bins-1]){fprintf(stderr,"Error: bin_max < r[num_bins-1].\n");exit(2);}

  //Acquire delta sigma from the block
  status |= c_datablock_get_array_ndim(block,delta_sigma_name,"delta_sigma",&ndims);
  status |= c_datablock_get_double_array_shape(block,delta_sigma_name,"delta_sigma",ndims,extents);
  if(status){fprintf(stderr,"Error on reading in delta_sigma (2D) dimensions.\n");exit(status);}
  int numsigr0 = extents[0],numsigr1=extents[1];
  double delta_sigma[numsigr0][numsigr1];
  status |= c_datablock_get_double_array(block,delta_sigma_name,"delta_sigma",(double *)delta_sigma,ndims,extents);
  if(status){fprintf(stderr,"Error on reading in delta_sigma.\n");exit(status);}

  //Declare and allocate the resulting array for the average surface density
  int ndimsd=2;
  int extents_out[]={numz,num_bins};
  double**ave_delsig;
  ave_delsig = (double**)malloc(numz*sizeof(double*));
  for(i=0;i<numz;i++)
    ave_delsig[i] = (double*)malloc(num_bins*sizeof(double));
  
  //Do the surface overdensity calculation
  //Note: we have to flatten the delta_sigma array
  status |= calc_ave_delsig(z,numz,r,numr,bin_edges,num_bins,(double *)delta_sigma,ave_delsig);

  //We have to flatten the ave_delsig array in order to write to the block
  double ave_delsig_out[numz][num_bins];
  for(i=0;i<numz;i++)
      for(j=0;j<num_bins;j++)
	ave_delsig_out[i][j]=ave_delsig[i][j];
  
  //Write ave_delsig to the block
  char *name = "ave_delta_sigma";
  status |= c_datablock_put_double_array(block,name,name,(double *)ave_delsig_out,ndimsd,extents_out);
  if(status){fprintf(stderr,"Error on writing out ave_delta_sigma to the block.\n");exit(status);}

  //Write the bin edges to the block
  int ndimedges=1;
  int extents_edges[]={num_bins+1};
  status |= c_datablock_put_double_array(block,name,"bin_edges",(double *)bin_edges,ndimedges,extents_edges);
  if(status){fprintf(stderr,"Error on writing out ave_delta_sigma to the block.\n");exit(status);}


  //Free the ave_delsig array
  for(i=0;i<numz;i++)
    free(ave_delsig[i]);
  free(ave_delsig);

  return 0;//No likelihood from this module
}

//In cleanup we have to free the configuration data.
int cleanup(void * config_in){
  ave_delsig_config * config = (ave_delsig_config *)config_in;
  free(config);
  return 0;
}
