#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int calc_sigma_r(double M,double r_scale,double delta,
		 double *z,int numz,
		 double *r,int numr,
		 double *xi_in, double omega_m,double H0, 
		 double **sigma_r,double offset);
