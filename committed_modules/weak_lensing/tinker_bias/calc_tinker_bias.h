#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int calc_tinker_bias(double delta,double M,
		     double*z,int numz,double*k,int numk,
		     double*PK_in,double H0,double omega_m,
		     double*bias,double*sigma2m);
