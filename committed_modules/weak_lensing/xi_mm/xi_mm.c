/*
  This module will calculate the mm correlation function.
 */
#include "cosmosis/datablock/c_datablock.h"
#include "cosmosis/datablock/section_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "calc_xi_mm.h"

//These are some section names we will use
const char * cosmo = COSMOLOGICAL_PARAMETERS_SECTION;
const char * dist = DISTANCES_SECTION;
const char * setup_name = "setup_delta_sigma";//Hard coded, for now

/* This struct will hold the data from the options.*/
typedef struct xi_mm_config{
  char *matter_power;//The kind of power spectrum
} xi_mm_config;

/* In setup, the options are obtained and passed to execute via the config object.*/
void * setup(c_datablock * options){
  DATABLOCK_STATUS status=0;
  xi_mm_config * config = malloc(sizeof(xi_mm_config));
  status |= c_datablock_get_string(options,OPTION_SECTION,"matter_power",&(config->matter_power));

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
  xi_mm_config * config = (xi_mm_config *)config_in;
  char *matter_power_blockname = config->matter_power;

  //Acquire the radii to evaluate the correlation function at.
  double *r;
  int numr;
  status |= c_datablock_get_double_array_1d(block,setup_name,"radii",&r,&numr);

  //for(i=0;i<numz;i++)
  //  printf("z=%e\n",z[i]);
  //double sigma8,sigma8_in;
  //status |= c_datablock_get_double(block,cosmo,"sigma_8",&sigma8);
  //status |= c_datablock_get_double(block,cosmo,"sigma8_input",&sigma8_in);
  //printf("\nsigma8    = %e\n",sigma8);
  //printf("sigma8_in = %e\n\n",sigma8_in);

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
  status |= c_datablock_get_double_array(block,matter_power_blockname,"p_k",(double *)PK,ndims,extents);
  if(status){fprintf(stderr,"Error on reading in PK.\n");exit(status);}
  if(numk!=numPK2||numz!=numPK1){fprintf(stderr,"Error, numz!=numPk1 or numk!=numPK2.");exit(1);}//Check sizes of arrays

  //Declare and allocate the final array to be written to the block
  int ndimxi2h = 2;
  int extentsxi2h[] = {numz,numr};
  double **xi_mm;
  xi_mm = (double **)malloc(numz*sizeof(double*));
  for(i=0;i<numz;i++){
    xi_mm[i]= (double *)malloc(numr*sizeof(double));
  }

  //Do the xi mm calculation
  //Note: we have to flatten the PK array
  status |= calc_xi_mm(r,numr,z,numz,k,numk,(double *)PK,numPK1,numPK2,xi_mm);
  
  //We have to flatten the xi_mm array in order to write to the block
  //This might not be necessary once the ND array stuff in CosmoSIS is implemented
  double xi_mm_out[numz][numr];
  for(i=0;i<numz;i++)
    for(l=0;l<numr;l++)
      xi_mm_out[i][l]=xi_mm[i][l];
  
  //Write xi_mm to the block
  char *name = "xi_mm";
  status |= c_datablock_put_double_array(block,name,name,(double *)xi_mm_out,ndimxi2h,extentsxi2h);
  if(status){fprintf(stderr,"Error on writing out xi_mm to the block.\n");exit(status);}

  //Free the xi_mm array
  for(i=0;i<numz;i++)
    free(xi_mm[i]);
  free(xi_mm);

  return 0;//No likelihood from this module
}

//In cleanup we have to free the configuration data.
int cleanup(void * config_in){
  xi_mm_config * config = (xi_mm_config *)config_in;
  free(config);
  return 0;
}
