/**
 * @file kde.c
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


#include "kde.h"

/** Estimate bandwidth using Silverman's "rule of thumb" 
 * (Silverman 1986, pg 48 eq 3.31).  This is the default
 * bandwith estimator for the R 'density' function.  */
double nrd0(double x[], const int N)
{
	gsl_sort(x, 1, N);
	double hi = gsl_stats_sd(x, 1, N);
	double iqr =
		gsl_stats_quantile_from_sorted_data (x,1, N,0.75) - 
        gsl_stats_quantile_from_sorted_data (x,1, N,0.25);
	double lo = GSL_MIN(hi, iqr/1.34);
	double bw = 0.9 * lo * pow(N,-0.2);
	return(bw);
}

/* kernels for kernel density estimates */
double gauss_kernel(double x)
{ 
	return exp(-(gsl_pow_2(x)/2))/(M_SQRT2*sqrt(M_PI)); 
}

double epanechnikov_kernel(double x)
{
  //fprintf(stderr, "x: %f\n", x);
  if( x > 1 || x < -1)
    {
      return 0.;
    }
  return 0.75 * (1. - x*x);
}

void  kerneldensity(double *samples, double *probs, double bandwith, double* nobs, size_t n, size_t npoints, int mode)
{
  size_t i,j;
  double h;
  h = GSL_MAX(nrd0(samples, n), 1e-6);
  fprintf(stderr, "bw: %f\n", h);
  if(bandwith > 0)
    h = bandwith;
  
  double prob = 0;
  
  
  for(j = 0; j < npoints; ++j){
    prob=0.;
    for(i=0; i < n; i++)
      {
	if( mode == 1) // gaussian
	  prob += gauss_kernel( (samples[i] - nobs[j])/h)/(n*h);
	if( mode == 2) // epanechnikov
	  prob += epanechnikov_kernel((samples[i] - nobs[j])/h/sqrt(5))/(n*h*sqrt(5));
      }
    probs[j] = prob;
  }
  
}

