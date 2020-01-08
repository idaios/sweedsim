#include "kde.h"


int main(int argc, char **argv){
  const gsl_rng_type * T;
  gsl_rng * r;
  int i, n = 1000, denPoints = 512;
  double a = -1, b = 1;
  double bw = -1;
  double *v = calloc(n, sizeof(double)), *density = calloc(denPoints+1, sizeof(double));
  
  /* create a generator chosen by the
     environment variable GSL_RNG_TYPE */
  
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  FILE *data = fopen("data.txt",  "w");
  
  /* print n random variates chosen from
     the poisson distribution with mean
     parameter mu */
  //hat_linear(v, 1000, density, 1, a, b, denPoints);

  for (i = 0; i < n; i++)
    {
      v[i] = gsl_ran_flat (r, a, b);
      fprintf (data, "%f\n", v[i]);
    }
  fclose(data);

  double* xpoints = calloc(denPoints+1, sizeof(double));
  xpoints[0] = a;
  double step = (double)(b - a)/(denPoints);
  for(i = 1; i < denPoints+1; ++i){
    xpoints[i] = xpoints[i-1] + step;
    //printf("point %d %f: ", i,xpoints[i]);
  }
  //printf("\n");

  kerneldensity(v, density, bw, xpoints, 1000, denPoints+1, 2);


  double sum = 0.;

  for (i = 0; i < denPoints; i++){
    printf ("%f\t%f\n", xpoints[i], density[i]);
    sum += density[i];
  }
  //printf("\ndensity: %f\n", sum);
  
  return 0;
  
}


