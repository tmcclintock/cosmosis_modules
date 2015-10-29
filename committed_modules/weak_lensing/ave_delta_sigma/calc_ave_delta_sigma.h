#include "gsl/gsl_spline.h"
#include "gsl/gsl_integration.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

int calc_ave_delsig(double *z,int numz,
		    double *r,int numr,
		    double *bin_edges,int num_bins,
		    double *delta_sigma_in,
		    double **ave_delsig);
