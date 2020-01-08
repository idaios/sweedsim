/**
 * @file kde.h
 * @author Carl Boettiger, <cboettig@gmail.com>
 * @section DESCRIPTION
 * Estimates the kernel density p(x) at a given value x from
 * an array of sample points.  Uses the default algorithm from
 * the R langauge's 'density' function.  Requires the GSL statistics
 * library.  
 *   
 * @section LICENCE
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */



#include <stdio.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>

void  kerneldensity(double *samples, double *probs, double bandwith, double* nobs, size_t n, size_t npoints, int mode);
