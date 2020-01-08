
/*
  
  The program has been written by Dick Hudson to generate the forward-in-time trajectory of the beneficial mutation. 
  A sister program called stepfn and stepfn2 is called to generate the trajectory backwards in time. 

  The program has been modified by Pavlos Pavlidis in 2015 to accept further arguments for the demographic model.   

  Below are the notes from Dick Hudson:
  -------------------------------------
  
  This program generates frequency trajectories of an allele that has freqeuency currently in the
  interval { pf-eps, pf+eps}.   It assumes that the mutation arose once, at a random time in the past, and
  no further mutation occurs.
  The command line arguments are nreps  genmax  s  h pf eps Npres  and seed.
  nreps:  number of trajectories to generate
  genmax:  It is assumed that the favored mutation occurred in the interval (0, genmax).
  s:      selection coefficient (  fitnesses:  1, 1+sh, 1+s ). 
  h:      dominance coefficient
  npres:     current diploid population size
  pf:     the final frequency attained by the favored allele (+/- eps). 
  seed:   for the random number generator.
  
  The function popsize(j), returns the diploid population size at time j generations in the past.  The user
  must supply this function.  Examples are provided.
  
  A Wright-Fisher model is assumed.  
  The first line of the output is nreps.
  The following line contains:   npoints:  n1   ( where n1 is the number of time points until present.)
  Following that are n1 lines where each line contains two numbers, the time point(= generation/(4*N) ), and
  the frequency.
  Time is measured in units of 4N generations.
  
  Compilation:  gcc -o trajdemog trajdemog.c popsize.c  binomial.c  rand1.c -lm
  usage:   trajdemog  100  2000  .02 0.5 0.8 .005 20000  931 >my.out
  output:
  
  trajdemog 100 2000 .02 0.5 0.8 .005 200000 931
  100
  #
  0.000000	0.000005
  0.000003	0.000005
  0.000005	0.000005
  0.000008	0.000010
  0.000010	0.000010
  0.000013	0.000020
  0.000015	0.000015
  ...
*/



/* 
   Modifications by Pavlos Pavlidis; November 2015
   
   -- The code accepts now the -eN flag to specify a stepwise demographic model
   -- All population sizes and generations are integers
   -- Support the -I flag; however the beneficial mutation trajectory is generated only
   in one population. 
   -- Makefile
   -- github repository
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#define MAXCH 100
#define TOL 1e-10
#define MAXXX 1000000

// the present day population size
extern int npres ;
extern double r ;



int main(int argc, char **argv )
{
  
  
  /* default values for the selection coefficient s,
     dominant coefficient h,
     current frequency of the beneficial mutation p,
     final frequency of the beneficial mutation pfinal,
     average fitness wbar,
     trajectory of the beneficial mutation freq
     tolerance level eps
  */
  double *s=NULL, h=0.5, p = 0., *pfinal=NULL, wbar, pprime, **freq, *eps=NULL;
  int nreps, gen, i , genmax=0, ngen, popmax, ch = 0, sweeptime = 0, maincounter = 0, poporigin = 0  ;
  float bnldev(float , int);
  int  nsuccess;
  int popsize(int );
  /* this function sorts a vector b based on a */
  int* intSortBonAint(int *a, int *b, int n);  
  int popsizeEN(int *sortedPopChanges, int *sortedTimeChanges, int n, int t, int npres);
  double ran1() ;
  unsigned short seedv[3], *seed48()  ;
  
  int *unsortedTimeChanges = calloc(MAXCH, sizeof(int));
  int *unsortedPopChanges = calloc(MAXCH, sizeof(int));
  int *timeChanges = calloc(MAXCH, sizeof(int));
  int *popChanges = calloc(MAXCH, sizeof(int));
  double *migvec = NULL;
  int npop = 0;

  if( argc < 6 ){
    fprintf(stderr,"trajdemog  -nreps -genmax -s -h  -pfinal -eps -npres -seed -eN\n");
    exit(1);
  }
  
   printf("// ");
  for( i=0; i<argc ; i++) printf(" %s",argv[i]);
  printf("\n");

  seedv[0] = 3579 ;

  for( i = 1; i < argc; ++i)
    {
      if(strcmp("-t", argv[i] ) == 0)
	{
	  assert(npop > 0);
	  poporigin = atoi(argv[++i]);
	  assert(poporigin < npop);
	  sweeptime = strtol(argv[++i], NULL, 10);
	  assert(sweeptime > 0);
	  continue;
	}

      if( strcmp("-nreps", argv[i]) == 0)
	{
	  nreps = strtol(argv[++i],NULL,10);
	  continue;
	}
      if( strcmp("-genmax", argv[i]) == 0)
	{
	  genmax = strtol(argv[++i], NULL, 10);
	  continue;
	}

      if( strcmp("-s", argv[i]) == 0 )
	{
	  assert(npop > 0);
	  s = calloc(npop, sizeof(double));
	  
	  for(int j = 0; j < npop; ++j)
	    s[j] = strtod( argv[++i], NULL);
	  
	  continue;
	}
      
      if( strcmp("-h", argv[i]) == 0 )
	{
	  h = strtod( argv[++i], NULL);
	  continue;
	}

      if(strcmp( "-pfinal", argv[i]) == 0)
	{
	  assert(npop > 0);
	  pfinal = calloc(npop, sizeof(double));
	  for( int j = 0; j < npop; ++j)
	    pfinal[j] = strtod( argv[++i], NULL);
	  continue;
	}

      if(strcmp("-eps", argv[i]) == 0 )
	{
	  assert(npop > 0);
	  eps = calloc(npop, sizeof(double));
	  for( int j = 0; j < npop; ++j)
	    eps[j] = strtod( argv[++i], NULL );
	  continue;
	}

      if(strcmp("-npres", argv[i]) == 0)
	{
	  npres = strtol( argv[++i], NULL, 10);
	  continue;
	}

      if(strcmp("-eN", argv[i]) == 0)
	{
	  unsortedTimeChanges[ch] = atoi(argv[++i]); //strtod(argv[++i], NULL);
	  unsortedPopChanges[ch] = atoi(argv[++i]);
	  ++ch;
	  continue;
	}

      if(strcmp("-mig", argv[i]) == 0)
	{
	  assert(npop > 0);
	  migvec = calloc(npop, sizeof(double));
	  for(int j = 0; j < npop; ++j)
	    migvec[j] = atof(argv[++i]);
	  continue;
	}
      
      if(strcmp("-npop", argv[i]) == 0)
	{
	  npop = atoi(argv[++i]);
	  continue;
	}
      
      if(strcmp("-seed", argv[i]) == 0)
	{
	  seedv[0] = strtol( argv[++i], NULL, 10);
	  continue;
	}
      
      fprintf(stderr, "Argument %s does not exist\n", argv[i]);
      exit(0);
    }

  assert(npop > 0);
  assert( sweeptime == 0 || genmax == 0 );
  timeChanges = intSortBonAint(unsortedTimeChanges, unsortedTimeChanges, ch);
  popChanges = intSortBonAint( unsortedTimeChanges, unsortedPopChanges, ch);
 
  seedv[1] = 27011;
  seedv[2] = 59243 ;
  seed48( seedv);


  popmax = 0 ;
  
  for( i=0; i<genmax ; i++) 
    if( popsizeEN(popChanges, timeChanges, ch, i, npres) > popmax ) 
      {
	popmax = popsizeEN(popChanges, timeChanges, ch, i, npres) ;
      }

  /* fprintf(stderr, "popmax: %d\n", popmax); */
  
  printf("%d\n",nreps); 
  
  nsuccess = 0 ;

  freq = (double**)calloc(npop, sizeof(double*));
    
  if(sweeptime > 0)
    {
      ngen = sweeptime;
      for(i = 0; i < npop; ++i)
	freq[i] = (double *)malloc( (unsigned)(ngen+1)*sizeof(double) ) ;
    }
  else
    {
      ngen = ran1()*genmax + 1 ;
      poporigin = 0;
      for(i = 0; i < npop; ++i)
	freq[i] = (double *)malloc( (unsigned)(genmax+1)*sizeof(double) ) ;
    }  

  
  
  /* char **allstrings = calloc( ngen, sizeof(char *) ); */

  /* for( i = 0; i < ngen; ++i) */
  /*   allstrings[i] = calloc(100, sizeof(char)); */


  //fprintf(stderr, "ngen: %d\n", ngen);
  
  
  while( nsuccess < nreps ){
    
    if(sweeptime > 0)
      ngen = sweeptime;
    else
      ngen = ran1()*genmax + 1 ;
    
    if( sweeptime > 0 || ran1() < popsizeEN(popChanges, timeChanges, ch, ngen, npres)/(double)popmax ) {

      for(int j = 0; j < npop; ++j)
	for(i = 0; i < ngen; ++i)
	  freq[j][i] = 0;
      /* gen = 0 ; */
      /* freq[0] = 0.0 ; */
      gen = 1 ;
      // popsizeEN returns the current effective population size
      freq[poporigin][1] = 5.0/(2.0*popsizeEN(popChanges, timeChanges, ch, ngen, npres) ) ;
      assert( freq[poporigin][1] > 0);
      assert( freq[poporigin][1] <= 1);

      // this loop runs over the generations
      while( (gen < ngen)  ){ //&& ( p > 0) && ( p<1.0) ){

	/* if(gen == 100) */
	/*   exit(1); */
	
	++gen;

	int condition = 1;
	for(int j = 0; j < npop; ++j){
	  p = 0.;
	  for( int jj = 0; jj < npop; ++jj)
	    if( jj != j)
	      p += (1.-migvec[j])*freq[j][gen-1] // the probability of no migration
		+ (migvec[j]/(npop-1)) * freq[jj][gen-1]; // the probability of migrants
	  
	  //fprintf(stderr, "pop: %d, gen: %d,  p: %.12f\n", j, gen, p); 
	  wbar = p*p*(1.+s[j]) + 2.*p*(1.-p)*(1.+s[j]*h)  + (1.-p)*(1.-p) ;
	  pprime = ( p*p*(1.+s[j]) + p*(1.-p)*(1.+s[j]*h) ) / wbar ;
	  //fprintf(stderr, "pop: %d, gen: %d,  p: %.12f, pprime:%.12f\n", j, gen, p, pprime);
	  double currentpop = popsizeEN(popChanges, timeChanges, ch, ngen-gen, npres);
	  
	  assert( pprime >= 0);
	  assert( pprime <= 1.0);
	  int successes = bnldev(pprime, 2*((int)currentpop));
	  
	  //fprintf(stderr, "successes: %d, N: %d\n", successes, 2*(int)currentpop);
	  /* fprintf(stderr, "pop: %d, pprime: %f, sel: %f, wbar: %f\n", j, pprime, s[j], wbar); */
	  
	  p = successes/2.0/currentpop;
	  /* fprintf(stderr, "pfinal: %.12f\n", p); */
	  
	  assert(gen < MAXXX);
	  freq[j][gen] = p ;

	  if(j == 0)
	    condition = condition && (p < 1 && p > 0);
	  /* fprintf(stderr, "condition: %d\n", condition); */
	}
	
	if(!condition)
	  break;
      }

      
      int conditionsMet = 0;
      for( int j = 0; j < npop; ++j)
	{
	  p = freq[j][gen];

	  if( (p >= pfinal[j] -eps[j]) && (p <= pfinal[j] + eps[j]) ){
	    conditionsMet++;
	  }
	}
      
      if(conditionsMet == npop)
	{
	  nsuccess++;
	  
	  if( (nreps > 9) && (nsuccess % (int)(0.1*nreps) == 0 ) ) 
	    fprintf(stderr,". ");
	  
	  printf("#\n");
      	
	  for( i=0; i<ngen; i++) 
	    {
	      printf("%6.12lf", i/(4.0*npres)) ;
	      for(int j = 0; j < npop; ++j){
		int fixationreached = 0 ; 
		if( freq[j][i] >= 1.00-TOL )
		  fixationreached = 1;
		if(fixationreached == 0)
		  printf("\t%6.12lf", freq[j][i]) ;
		else
		  printf("\t%6.12lf", 1.0) ;
	      }
	      printf("\n");
	    }
	  
	  maincounter = 0;
	}
      else
	maincounter++;
      
      if(maincounter > 0 && (maincounter % 100000) == 0)
	{
	  fprintf(stderr, "warning... difficult to reach fixation... %d tries\n", maincounter);
	}
    }
        
  }
  return 1;
}
  
