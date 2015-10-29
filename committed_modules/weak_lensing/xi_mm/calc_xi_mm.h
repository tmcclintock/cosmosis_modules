#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int calc_xi_mm(double *radii,int numr,double *z,int numz,double *k,int numk,double *PKin,int numPK1,int numPK2,double **xi_mm);
