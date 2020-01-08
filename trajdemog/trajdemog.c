
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
  double s=0.01, h=0.5, p = 0., pfinal=1.0, wbar, pprime, *freq, eps=0.001;
  
  int nreps, gen, i , genmax=0, ngen, popmax, ch = 0, sweeptime = 0, maincounter = 0  ;
  
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
	  sweeptime = strtol(argv[++i], NULL, 10);
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
	  s = strtod( argv[++i], NULL);
	  continue;
	}
      
      if( strcmp("-h", argv[i]) == 0 )
	{
	  h = strtod( argv[++i], NULL);
	  continue;
	}

      if(strcmp( "-pfinal", argv[i]) == 0)
	{
	  pfinal = strtod( argv[++i], NULL);
	  continue;
	}

      if(strcmp("-eps", argv[i]) == 0 )
	{
	  eps = strtod( argv[++i], NULL );
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

      if(strcmp("-seed", argv[i]) == 0)
	{
	  seedv[0] = strtol( argv[++i], NULL, 10);
	  continue;
	}
      
      fprintf(stderr, "Argument %s does not exist\n", argv[i]);
      exit(0);
    }

  assert( sweeptime == 0 || genmax == 0 );

  timeChanges = intSortBonAint(unsortedTimeChanges, unsortedTimeChanges, ch);

  popChanges = intSortBonAint( unsortedTimeChanges, unsortedPopChanges, ch);
 
  /* nreps = strtol(argv[1],NULL,10); */

  /* genmax = strtol(argv[2],NULL,10); */

  /* s = strtod( argv[3],NULL ) ; */

  /* h = strtod( argv[4],NULL ) ; */

  /* pfinal = strtod( argv[5],NULL ) ; */

  /* eps = strtod( argv[6],NULL ) ; */
  
  /* npres = strtol(argv[7],NULL,10); */

  /* fprintf(stderr,"%d %lf %lf %lf %d\n",nreps, s, h,  pfinal, npres); */
  
      //if( argc > 8 ) seedv[0] = strtol( argv[8] ,NULL,10) ;
  
  /* fprintf(stderr,"seed0: %d\n",seedv[0] ) ; */
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

    
  if(sweeptime > 0)
    {
      ngen = sweeptime;
      freq = (double *)malloc( (unsigned)(ngen+1)*sizeof(double) ) ;
    }
  else
    {
      ngen = ran1()*genmax + 1 ;
      freq = (double *)malloc( (unsigned)(genmax+1)*sizeof(double) ) ;
    }
  

  
  
  /* char **allstrings = calloc( ngen, sizeof(char *) ); */

  /* for( i = 0; i < ngen; ++i) */
  /*   allstrings[i] = calloc(100, sizeof(char)); */

 
  
  while( nsuccess < nreps ){
    
    
    if(sweeptime > 0)
      ngen = sweeptime;
    else
      ngen = ran1()*genmax + 1 ;
    
    if( sweeptime > 0 || ran1() < popsizeEN(popChanges, timeChanges, ch, ngen, npres)/(double)popmax ) {

      for(i = 0; i < ngen; ++i)
	{
	  freq[i] = 0;
	  
	}

      gen = 0 ;
      freq[0] = 0.0 ;
      gen = 1 ;

      // popsizeEN returns the current effective population size
      p = freq[1] = 1.0/(2.0*popsizeEN(popChanges, timeChanges, ch, ngen, npres) ) ;
      assert( p >= 0);
      assert( p <= 1);

      // this loop runs over the generations
      while( (gen < ngen)  && ( p > 0) && ( p<1.0) ){
	gen++;
	wbar = p*p*(1.+s) + 2.*p*(1.-p)*(1.+s*h)  + (1.-p)*(1.-p) ;
	pprime = ( p*p*(1.+s) + p*(1.-p)*(1.+s*h) ) / wbar ; // next frequency
	
	double currentpop = popsizeEN(popChanges, timeChanges, ch, ngen-gen, npres);

	//assert(currentpop == npres);

	assert( pprime >= 0);

	assert( pprime <= 1.0);
	
	int successes = bnldev(pprime, 2*((int)currentpop));
	
	p = successes/2.0/currentpop;

	
	assert(gen < MAXXX);
	
	freq[gen] = p ;

      }
      
      if( (p > pfinal-eps) && (p<pfinal+eps) ){
	
	nsuccess++;
	
	if( (nreps > 9) && (nsuccess % (int)(0.1*nreps) == 0 ) ) 
	  fprintf(stderr,". ");
	
	printf("#\n");

	int fixationreached = 0;
	
	for( i=0; i<ngen; i++) 
	  {
	    
	    if(fixationreached == 0)
	      printf("%6.12lf\t%6.12lf\n", i/(4.0*npres), freq[i]) ;
	    else
	      printf("%6.12lf\t%6.12lf\n", i/(4.0*npres), 1.0) ;

	    if( freq[i] >= 1.00-TOL )
	      fixationreached = 1;
	      	    
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
  
  fprintf(stderr," !\n");

  return 1;
  
}
  
