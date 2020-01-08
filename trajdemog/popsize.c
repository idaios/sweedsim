/*
 *  popsize.c
 *  
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int na = 10000 ;
int nb = 1000 ;
int nc = 5000 ;
int t1 = 500 ;
int t2 = 2500 ;
int t3 = 200 ;
int npres = 100000 ;


double* sortBonA(double *a, double *b, int n)
{
  
  int c, d;

  double *k = calloc(n, sizeof(double)), swap;

  for( c = 0; c < n; ++c)
    k[c] = b[c];
  
  for (c = 0 ; c < ( n - 1 ); c++)
    {
      for (d = 0 ; d < n - c - 1; d++)
	{
	  if (a[d] > a[d+1]) 
	    {
	      swap       = k[d];
	      k[d]   = k[d+1];
	      k[d+1] = swap;
	    }
	}
    }
  
  return k;
}


int* sortBonAint(double *a, int *b, int n)
{
  
  int c, d;

  int *k = calloc(n, sizeof(int)), swap;

  for( c = 0; c < n; ++c)
    k[c] = b[c];
  
  for (c = 0 ; c < ( n - 1 ); c++)
    {
      for (d = 0 ; d < n - c - 1; d++)
	{
	  if (a[d] > a[d+1]) 
	    {
	      swap       = k[d];
	      k[d]   = k[d+1];
	      k[d+1] = swap;
	    }
	}
    }
  
  return k;
}

int* intSortBonAint(int *a, int *b, int n)
{
  
  int c, d;

  int *k = calloc(n, sizeof(int)), swap;

  for( c = 0; c < n; ++c)
    k[c] = b[c];
  
  for (c = 0 ; c < ( n - 1 ); c++)
    {
      for (d = 0 ; d < n - c - 1; d++)
	{
	  if (a[d] > a[d+1]) 
	    {
	      swap   = k[d];
	      k[d]   = k[d+1];
	      k[d+1] = swap;
	    }
	}
    }
  
  return k;
}



int popsizeEN(int *sortedPopChanges, int *sortedTimeChanges, int n, int t, int npres)
{
  int i, pop = 999999;

  if( n == 0)
    return npres;

  for( i = 0; i < n; ++i)
    {
      
      if( t < sortedTimeChanges[i] )
	{
	  
	  if( i == 0 )
	    {
	      pop =  npres;
	      break;
	    }
	  
	  
	  pop = sortedPopChanges[i-1];
	  
	  break;
	  
	}	  
    }

  if( i == n)
    pop = sortedPopChanges[i-1];
  
  return pop;
  
}

	int
popsize(int t)
{
double r ;
int pop ;

	r = (1.0/t1)*log(npres/(double)nc) ;
	if( t <= t1 ) pop = (int) ( npres*exp( -r*t ) ) ;
	else if (t <= t2+t1 ) pop = nc ;
	else if ( t <= t1+t2+t3 ) pop = nb ;
	else pop = na ;
	return(pop);
}

