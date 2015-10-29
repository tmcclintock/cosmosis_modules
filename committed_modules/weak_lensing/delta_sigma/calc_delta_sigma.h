#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int calc_delta_sigma(double M,double r_scale,double delta,
		     double omega_m,double H0,
		     double *z,int numz,double *r,int numr,
		     double *sigma_r_in,double **delta_sigma);
