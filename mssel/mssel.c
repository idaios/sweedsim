/***** mssel.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
*	Usage is shown by typing mssel without arguments.   
        usage: mssel nanc nder howmany selfilename selspot  -t  theta   -r 4Nr  nsites  [options]

	   nanc is the number of gametes per sample with the ancestral allele at the selected site.
	   nder is the number of gametes per sample with the derived allele at the selected site.
	   howmany is the number of samples to produce.
	   selfilename is the file with the allele frequencies through time in each subpopulation.
	   selspot is the position of the selected site.

 selfile format:
nreps
maxpop  
n: nfreqtimepoints
t0  f1 f2 ... fmaxpop
t1  f1 f2 ... fmaxpop
...
tn-1 f1 f2 ... fmaxpop
n: nfreqtimepoints
t0  f1 f2 ... fmaxpop
t1  f1 f2 ... fmaxpop
...
tn-1 f1 f2 ... fmaxpop
...
n: nfreqtimepoints
t0  f1 f2 ... fmaxpop
t1  f1 f2 ... fmaxpop
...
tn-1 f1 f2 ... fmaxpop

	To compile:  cc -o mssel  mssel.c  streecsel.c  rand1.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
		some other optimization switches might be usefully employed with some 
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).


***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "mssel.h"
#include "kde.h"

#define SITESINC 10 

static unsigned maxsites = SITESINC ;

struct segl {
  int beg;
  struct _mssel_node *ptree;
  int next;
};

double *posit ;
double segfac ;
int count, ntbs, nseeds  ;
struct _mssel_params pars ;	

#define STATE1 '1'
#define STATE2 '0'

int _mssel_pickb(int nsam, struct _mssel_node *ptree, double tt)
{
  double x, y, _mssel_ran1();
  int i;

  x = _mssel_ran1()*tt;
  for( i=0, y=0; i < 2*nsam-2 ; i++) {
    y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
    if( y >= x ) return( i ) ;
  }
  return( i );
}



int _mssel_pickbmf(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct _mssel_node *ptree;
	double tt;
{
	double x, y, _mssel_ran1();
	int i;

	x = _mssel_ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	  if( y >= x ) return( i ) ;
	}
	return( i );
}




int _mssel_cmpfunc (const void * a, const void * b) {
  if( *(double*)a - *(double*)b > 0)
    return 1;
  return -1;
}




/* a slight modification of crecipes version */
double _mssel_gasdev(double m, double v)
{
  int iset=0;
  static float gset;
  float fac,r,v1,v2;
  double _mssel_ran1();
  //fprintf(stderr, "%d\n", iset);
  if  (iset == 0) {
    do {
      v1=2.0*_mssel_ran1()-1.0;
      v2=2.0*_mssel_ran1()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);
    gset= v1*fac;
    iset=1;
    return( m + sqrt(v)*v2*fac);
  } else {
    iset=0;
    return( m + sqrt(v)*gset ) ;
  }
}



void _mssel_mnmial(n,nclass,p,rv)
	int n, nclass, rv[];
	double p[];
{
  	double _mssel_ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = _mssel_ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return;
}

void _mssel_parens( struct _mssel_node *ptree, int *descl, int *descr,  int noden)
{
	double time ;

   if( descl[noden] == -1 ) {
	printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	}
   else{
	printf("(");
	_mssel_parens( ptree, descl,descr, descl[noden] ) ;
	printf(",");
	_mssel_parens(ptree, descl, descr, descr[noden] ) ;
	if( (ptree+noden)->abv == 0 ) printf(");\n"); 
	else {
	  time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
	  printf("):%5.3lf", time );
	  }
        }
}


void _mssel_prtree( ptree, nsam)
	struct _mssel_node *ptree;
	int nsam;
{
	int i, *descl, *descr ;
	//void _mssel_parens( struct _mssel_node *ptree, int *descl, int *descr, int noden );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	_mssel_parens( ptree, descl, descr, 2*nsam-2);
	free( descl ) ;
	free( descr ) ;
}


double _mssel_ttimemf( ptree, nsam, mfreq)
	struct _mssel_node *ptree;
	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}




/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

double _mssel_ttime( ptree, nsam)
	struct _mssel_node *ptree;
	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
	return(t);
}


int _mssel_tdesn(ptree, tip, node )
     struct _mssel_node *ptree;
	int tip, node;
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}


void _mssel_make_gametes(int nsam, int mfreq, struct _mssel_node *ptree, double tt, int newsites, int ns, char **list )
{
  int  tip, j,  node ;
  //int pickb(int nsam, struct _mssel_node *ptree, double tt);
  
  //    pickbmf(int nsam, int mfreq, struct _mssel_node *ptree, double tt) ;
  
  for(  j=ns; j< ns+newsites ;  j++ ) {
    if( mfreq == 1 ) node = _mssel_pickb(  nsam, ptree, tt);
    else node = _mssel_pickbmf(  nsam, mfreq, ptree, tt);
    for( tip=0; tip < nsam ; tip++) {
      if( _mssel_tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
      else list[tip][j] = STATE2 ;
    }
  }
}


int _mssel_poisso(double u)
{
  
  
  double  cump, ru, _mssel_ran1(), p; //, gasdev() ;
  int i=1;
  int tmp = (int)(0.5 + _mssel_gasdev(u,u));
  if( u > 30. )
    return tmp;
  ru = _mssel_ran1();
  p = exp(-u);
  if( ru < p) return 0;
  cump = p;
  while( ru > ( cump += (p *= u/i ) ) )
    i++;
  return(i);
}

void _mssel_ranvec(int n,double pbuf[])
{
	int i;
	double _mssel_ran1();
	
	for(i=0; i<n; i++)
	  pbuf[i] = _mssel_ran1();
	
	return;
}

void _mssel_order(int n, double pbuf[])
{
  int gap, i, j;
  double temp;
  
  for( gap= n/2; gap>0; gap /= 2)
    for( i=gap; i<n; i++)
      for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
	temp = pbuf[j];
	pbuf[j] = pbuf[j+gap];
	pbuf[j+gap] = temp;
      }
}



void _mssel_ordran(int n, double pbuf[])
{
  _mssel_ranvec(n,pbuf);
  _mssel_order(n,pbuf);
  return;
}

void _mssel_locate(n,beg,len,ptr)
	int n;
	double beg, len, *ptr;
{
  int i;
  
  _mssel_ordran(n,ptr);
  for(i=0; i<n; i++)
    ptr[i] = beg + ptr[i]*len ;
}


void  _mssel_ndes_setup(struct _mssel_node *ptree, int nsam )
{
  int i ;
  for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
  for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
  for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;
}

void _mssel_biggerlist(int nsam,  char** list )
{
  int i;
  
  /*  fprintf(stderr,"maxsites: %d\n",maxsites);  */	
  for( i=0; i<nsam; i++){
    list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
    if( list[i] == NULL ) perror( "realloc error. bigger");
  }
}

int _mssel_gensam( char **list, double *pprobss, double *ptmrca, double *pttot , int *pselsite ) 
{


  int nsegs, i, k, j, seg, ns, start, end, len, segsit , ind, pop, allele, selsite ;
  struct segl *seglst, *_mssel_segtre_mig(struct c_params *p, int *nsegs ) ; /* used to be: [MAXSEG];  */
  double nsinv,  tseg, tt; //, ttime(struct _mssel_node *, int nsam), ttimemf(struct _mssel_node *, int nsam, int mfreq) ;
  double *pk;
  int *ss;
  int segsitesin,nsites;
  double theta, es, selpos  ;
  int nsam, mfreq ;
  //	void prtree( struct _mssel_node *ptree, int nsam);
  //void make_gametes(int nsam, int mfreq,  struct _mssel_node *ptree, double tt, int newsites, int ns, char **list );
  //void _mssel_ndes_setup( struct _mssel_node *, int nsam );
  nsites = pars.cp.nsites ;
  nsinv = 1./nsites;
  seglst = _mssel_segtre_mig(&(pars.cp),  &nsegs ) ;
  
  nsam = pars.cp.nsam;
  segsitesin = pars.mp.segsitesin ;
  theta = pars.mp.theta ;
  mfreq = pars.mp.mfreq ;
  
  if( pars.mp.treeflag ) {
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
	end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
	start = seglst[seg].beg ;
	len = end - start + 1 ;
	fprintf(stdout,"[%d]", len);
      }
      _mssel_prtree( seglst[seg].ptree, nsam ) ;
      if( (segsitesin == 0) && ( theta == 0.0 ) && ( pars.mp.timeflag == 0 ) ) 
	free(seglst[seg].ptree) ;
    }
  }
  
  if( pars.mp.timeflag ) {
    tt = 0.0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      if( mfreq > 1 ) _mssel_ndes_setup( seglst[seg].ptree, nsam );
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      if( (nsegs==1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) )
	*ptmrca = (seglst[seg].ptree + 2*nsam-2) -> time ;
      len = end - start + 1 ;
      tseg = len/(double)nsites ;
      if( mfreq == 1 ) tt += _mssel_ttime(seglst[seg].ptree,nsam)*tseg ;
      else tt += _mssel_ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
      if( (segsitesin == 0) && ( theta == 0.0 )  ) 
	free(seglst[seg].ptree) ;
    }
    *pttot = tt ;
  }	
  
  if( (segsitesin == 0) && ( theta > 0.0)   ) {
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      if( mfreq > 1 ) _mssel_ndes_setup( seglst[seg].ptree, nsam );
      
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      len = end - start + 1 ;
      tseg = len*(theta/nsites) ;
      
      if( mfreq == 1) tt = _mssel_ttime(seglst[seg].ptree, nsam);
      else tt = _mssel_ttimemf(seglst[seg].ptree, nsam, mfreq );
      segsit = _mssel_poisso( tseg*tt );
      if( (segsit + ns +1) >= maxsites ) {
	maxsites = segsit + ns + SITESINC ;
	posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
	_mssel_biggerlist(nsam, list) ; 
      }

      _mssel_make_gametes(nsam,mfreq,seglst[seg].ptree,tt, segsit, ns, list );
      
      free(seglst[seg].ptree) ;
      _mssel_locate(segsit,start*nsinv, len*nsinv,posit+ns);   
      ns += segsit;
    }
  }
  else if( segsitesin > 0 ) {
    
    pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
    ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
    if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");
    
    
    tt = 0.0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      if( mfreq > 1 )
	_mssel_ndes_setup( seglst[seg].ptree, nsam );
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      len = end - start + 1 ;
      tseg = len/(double)nsites ;
      if( mfreq == 1 ) pk[k] = _mssel_ttime(seglst[seg].ptree,nsam)*tseg ;
      else pk[k] = _mssel_ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
      tt += pk[k] ;
    }
    if( theta > 0.0 ) { 
      es = theta * tt ;
      *pprobss = exp( -es )*pow( es, (double) (segsitesin-pars.selsitepolyflag) ) / segfac ;
    }
    if( tt > 0.0 ) {
      for (k=0;k<nsegs;k++) pk[k] /= tt ;
      _mssel_mnmial(segsitesin-pars.selsitepolyflag,nsegs,pk,ss);
    }
    else
      for( k=0; k<nsegs; k++) ss[k] = 0 ;
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      len = end - start + 1 ;
      tseg = len/(double)nsites;
      _mssel_make_gametes(nsam,mfreq,seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);
      
      free(seglst[seg].ptree) ;
      _mssel_locate(ss[k],start*nsinv, len*nsinv,posit+ns);   
      ns += ss[k] ;
    }
    
    free(pk);
    free(ss);
    
  }


  if( pars.selsitepolyflag){
    selsite = 0 ;
    selpos = (pars.cp.selspot+1)/(double)nsites ;
    while( (posit[selsite] < selpos) &&( selsite < ns)  ) selsite++;
    *pselsite = selsite ;
    for( i =ns-1; i >= selsite ; i-- ){
      posit[i+1] = posit[i] ;
      for( ind = 0; ind < nsam; ind++) list[ind][i+1] = list[ind][i];
    }
    posit[selsite] = selpos ;
    for(pop=ind=0;pop<pars.cp.npop;pop++){
      for(allele = 0; allele < 2 ; allele++){
	for(j=0; j<pars.cp.config[allele][pop];j++,ind++) {
	  list[ind][selsite] = '0' + allele  ;
	}
      }
    }
		ns++;
  }
  
  
  for(i=0;i<nsam;i++) list[i][ns] = '\0' ;
  return( ns ) ;
}


void _mssel_addtoelist( struct _mssel_devent *pt, struct _mssel_devent *elist ) 
{
  struct _mssel_devent *plast, *pevent, *ptemp  ;
  
  pevent = elist ;
  while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
    plast = pevent ;
    pevent = pevent->nextde ;
  }
  ptemp = plast->nextde ;
  plast->nextde = pt ;
  pt->nextde = ptemp ;
}


void  _mssel_free_eventlist( struct _mssel_devent *pt, int npop )
{
  struct _mssel_devent *next ;
  int pop ;
  
  while( pt != NULL){
    next = pt->nextde ;
    if( pt->detype == 'a' ) {
      for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
      free( pt->mat );
    }
    free(pt);
    pt = next ;
  }
}


void _mssel_usage()
{
  fprintf(stderr,"usage: mssel nsam howmany nanc nder selfilename selspot -r rho nsites [options] \n");
  fprintf(stderr,"  Options: \n"); 
  fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
  fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
  fprintf(stderr,"\t -T          (Output gene tree.)\n");
  fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
  fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
  fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
  fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
  fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
  fprintf(stderr,"\t -I npop n1anc n1der n2anc n2der ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
  fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
  fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
  fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
  fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
  fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
  fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
  fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
  fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
  fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
  fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
  fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
  fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
  fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
  fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
  fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
  fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
  fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
  fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");  
  fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");
  
  exit(1);
}



void _mssel_argcheck( int arg, int argc, char *argv[] )
{
  if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
    fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
    fprintf(stderr,"For usage type: mssel<return>\n");
    exit(0);
  }
}




static int NSEEDS = 3 ;

int  _mssel_getpars(int argc, char *argv[], int *phowmany )
{
  int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
  double migr, mij, psize, palpha ;
  //void addtoelist( struct _mssel_devent *pt, struct _mssel_devent *elist ); 
  //void _mssel_argcheck( int arg, int argc, char ** ) ;
  int _mssel_commandlineseed( char ** ) ;
  //void free_eventlist( struct _mssel_devent *pt, int npop );
  struct _mssel_devent *ptemp , *pt ;
  FILE *pf ;
  char ch3, line[1001] ;
  
  while(strcmp(argv[arg], "-mssel") != 0 && arg < argc - 1){
    arg++;	  
  }
  
  
  if( count == 0 ) {
    if(argc - arg < 6 ){
      return 1;
    }
    
    
    if( argc < 6 ){ fprintf(stderr,"Too few command line arguments\n"); _mssel_usage();}
    pars.cp.config = (int **)malloc( 2*sizeof( int *)) ;
    pars.cp.config[0] = (int *)malloc( sizeof(int)) ;
    pars.cp.config[1] = (int *)malloc( sizeof(int)) ;
    
    pars.cp.nsam = atoi( argv[arg+1] ) ;
    
    pars.cp.config[0][0] = atoi( argv[arg+3] );
    pars.cp.config[1][0] = atoi( argv[arg+4] );
    
    fprintf(stderr, "%d %d %d\n", pars.cp.nsam, pars.cp.config[0][0], pars.cp.config[1][0]);
    
    if( pars.cp.config[0][0] + pars.cp.config[1][0] != pars.cp.nsam ){
      fprintf(stderr,"Error.  nanc %d + nder %d is not equal to pars.cp.nsam %d.\n", pars.cp.config[0][0], pars.cp.config[1][0], pars.cp.nsam);
      assert(0);
      _mssel_usage();
      
    }
    
    if( pars.cp.nsam <= 0 ) { fprintf(stderr,"First argument error. nsam <= 0. \n"); _mssel_usage();}
    *phowmany = atoi( argv[arg+2] );
    if( *phowmany  <= 0 ) { fprintf(stderr,"Second argument error. howmany <= 0. \n"); _mssel_usage();}
    pars.selfile = fopen(argv[arg+5],"r");
    if( pars.selfile  == NULL ) {
      pars.selfile = stdin ;
      fprintf(stderr," trajectories from stdin\n" ); 
    }
    
    /* skip the lines starting with '#' */
    do {
      fgets( line, 1000, pars.selfile );
    }while ( line[0] == '#' );
    
    
    /* skip lines starting with '/'
       here we assume that some lines have been already read by the previous block of code
    */
    while( line[0] == '/' )
      {
	fgets( line, 1000, pars.selfile );
      }
    
    sscanf(line,"ntraj: %d",&(pars.cp.nrepsfreq)  );
    fscanf(pars.selfile,"npop: %d",&(pars.cp.maxpop) ); 
    fgets( line, 1000, pars.selfile );
    
    pars.cp.selspot = atoi( argv[arg+6] ) - 1 ;
    pars.commandlineseedflag = 0 ;
    pars.output_precision = 4 ;
    pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ;
    pars.cp.track_len = 0. ;
    pars.cp.npop = npop = 1 ;
    pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
    pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
    pars.cp.mig_mat[0][0] =  0.0 ;
    pars.mp.segsitesin = 0 ;
    pars.mp.treeflag = 0 ;
    pars.mp.timeflag = 0 ;
    pars.mp.mfreq = 1 ;
    pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
    (pars.cp.size)[0] = 1.0  ;
    pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
    (pars.cp.alphag)[0] = 0.0  ;
    pars.cp.nsites = 2 ;
  }
  else{
    pars.cp.config[0][0] = atoi( argv[arg+3] );
    pars.cp.config[1][0] = atoi( argv[arg+4] );
    pars.cp.selspot = atoi( argv[arg+6] ) - 1 ;  /* internally spots are 0 -> (nsites-1),  externally, 1 -> nsites */ 
    if( pars.cp.config[0][0] + pars.cp.config[1][0] != pars.cp.nsam ){
      fprintf(stderr,"nsam can't change.\n"); 
      _mssel_usage();
    }
    
    npop = pars.cp.npop ;
    _mssel_free_eventlist( pars.cp.deventlist, npop );
  }
  
  if( (pars.cp.config[0][0] > 0 ) && ( pars.cp.config[1][0] > 0 ) && (pars.cp.selspot >= 0) ) pars.selsitepolyflag = 1 ;
  else pars.selsitepolyflag = 0 ;
  
  pars.cp.deventlist = NULL ;
  
  arg = arg+7 ;
  
  while( arg < argc ){
    if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ? arg #: %d\n", argv[arg], arg+1 ); _mssel_usage();}
    switch ( argv[arg][1] ){
    case 'f' :
      if( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); exit(1); }
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pf = fopen( argv[arg], "r" ) ;
      if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", argv[arg] ); exit(0);}
      arg++;
      argc++ ;
      argv = (char **)malloc(  (unsigned)(argc+1)*sizeof( char *) ) ;
      argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
      argstart = arg ;
      while( fscanf(pf," %s", argv[arg]) != EOF ) {
	arg++;
	argc++;
	argv = (char **)realloc( argv, (unsigned)argc*sizeof( char*) ) ;
	argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
      }
      fclose(pf);
      argc--;
      arg = argstart ;
      break;
    case 'r' : 
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pars.cp.r = atof(  argv[arg++] );
      _mssel_argcheck( arg, argc, argv);
      pars.cp.nsites = atoi( argv[arg++]);
      if( pars.cp.nsites <2 ){
	fprintf(stderr,"with -r option must specify both rec_rate and nsites>1\n");
	_mssel_usage();
      }
      if( pars.cp.nsites < pars.cp.selspot ) pars.selsitepolyflag = 0 ;
      break;
    case 'p' :
      arg++;
      _mssel_argcheck(arg,argc,argv);
      pars.output_precision = atoi( argv[arg++] ) ;
      break;				
    case 'c' : 
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pars.cp.f = atof(  argv[arg++] );
      _mssel_argcheck( arg, argc, argv);
      pars.cp.track_len = atof( argv[arg++]);
      if( pars.cp.track_len <1. ){
	fprintf(stderr,"with -c option must specify both f and track_len>0\n");
	_mssel_usage();
      }
      break;		
    case 't' : 
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pars.mp.theta = atof(  argv[arg++] );
      break;
    case 's' : 
      arg++;
      _mssel_argcheck( arg, argc, argv);
      if( argv[arg-1][2] == 'e' ){  /* command line seeds */
	pars.commandlineseedflag = 1 ;
	if( count == 0 ) nseeds = _mssel_commandlineseed(argv+arg );
	arg += nseeds ;
      }
      else {
	pars.mp.segsitesin = atoi(  argv[arg++] );
      }
      break;
    case 'F' : 
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pars.mp.mfreq = atoi(  argv[arg++] );
      if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
	fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
	_mssel_usage();
      }
      break;
    case 'T' : 
      pars.mp.treeflag = 1 ;
      arg++;
      break;
    case 'L' : 
      pars.mp.timeflag = 1 ;
      arg++;
      break;
    case 'I' : 
      arg++;
      if( count == 0 ) {
	_mssel_argcheck( arg, argc, argv);
	pars.cp.npop = atoi( argv[arg]);
	pars.cp.config[0] = (int *) realloc( pars.cp.config[0], (unsigned)( pars.cp.npop*sizeof( int)));
	pars.cp.config[1] = (int *) realloc( pars.cp.config[1], (unsigned)( pars.cp.npop*sizeof( int)));
	npop = pars.cp.npop ;
      }
      arg++;
      for( i=0; i< pars.cp.npop; i++) {
	_mssel_argcheck( arg, argc, argv);
	pars.cp.config[0][i] = atoi( argv[arg++]);
	_mssel_argcheck( arg, argc, argv);
	pars.cp.config[1][i] = atoi( argv[arg++]);
      }
      if( count == 0 ){
	pars.cp.mig_mat = 
	  (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
	pars.cp.mig_mat[0] = 
	  (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
	for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] = 
					(double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
	pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
	pars.cp.alphag = 
	  (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
	for( i=1; i< pars.cp.npop ; i++) {
	  (pars.cp.size)[i] = (pars.cp.size)[0]  ;
	  (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
	}
      }
      if( (arg <argc) && ( argv[arg][0] != '-' ) ) {
	_mssel_argcheck( arg, argc, argv);
	migr = atof(  argv[arg++] );
      }
      else migr = 0.0 ;
      for( i=0; i<pars.cp.npop; i++) 
	for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
      for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
      break;
    case 'm' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); _mssel_usage();}
      if( argv[arg][2] == 'a' ) {
	arg++;
	for( pop = 0; pop <npop; pop++)
	  for( pop2 = 0; pop2 <npop; pop2++){
	    _mssel_argcheck( arg, argc, argv);
	    pars.cp.mig_mat[pop][pop2]= atof( argv[arg++] ) ;
	  }
	for( pop = 0; pop < npop; pop++) {
	  pars.cp.mig_mat[pop][pop] = 0.0 ;
	  for( pop2 = 0; pop2 < npop; pop2++){
	    if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
	  }
	}	
      }
      else {
	arg++;
	_mssel_argcheck( arg, argc, argv);
	i = atoi( argv[arg++] ) -1;
	_mssel_argcheck( arg, argc, argv);
	j = atoi( argv[arg++] ) -1;
	_mssel_argcheck( arg, argc, argv);
	mij = atof( argv[arg++] );
	pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j]  ;
	pars.cp.mig_mat[i][j] = mij;
      }
      break;
    case 'n' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); _mssel_usage();}
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pop = atoi( argv[arg++] ) -1;
      _mssel_argcheck( arg, argc, argv);
      psize = atof( argv[arg++] );
      pars.cp.size[pop] = psize ;
      break;
    case 'g' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); _mssel_usage();}
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pop = atoi( argv[arg++] ) -1;
      if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); _mssel_usage(); }
      palpha = atof( argv[arg++] );
      pars.cp.alphag[pop] = palpha ;
      break;
    case 'G' :
      arg++;
      if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); _mssel_usage(); }
      palpha = atof( argv[arg++] );
      for( i=0; i<pars.cp.npop; i++) 
	pars.cp.alphag[i] = palpha ;
      break;
    case 'e' :
      pt = (struct _mssel_devent *)malloc( sizeof( struct _mssel_devent) ) ;
      pt->detype = argv[arg][2] ;
      ch3 = argv[arg][3] ;
      arg++;
      _mssel_argcheck( arg, argc, argv);
      pt->time = atof( argv[arg++] ) ;
      pt->nextde = NULL ;
      if( pars.cp.deventlist == NULL ) 
	pars.cp.deventlist = pt ;
      else if ( pt->time < pars.cp.deventlist->time ) { 
	ptemp = pars.cp.deventlist ;
	pars.cp.deventlist = pt ;
	pt->nextde = ptemp ;	
      }	
      else
	_mssel_addtoelist( pt, pars.cp.deventlist ) ;
      switch( pt->detype ) {
      case 'N' :
	_mssel_argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'G' :
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); _mssel_usage(); }
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'M' :
	_mssel_argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'n' :
	_mssel_argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	_mssel_argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'g' :
	_mssel_argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); _mssel_usage(); }
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 's' :
	_mssel_argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	_mssel_argcheck( arg, argc, argv);
	pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'm' :
	if( ch3 == 'a' ) {
	  pt->detype = 'a' ;
	  _mssel_argcheck( arg, argc, argv);
	  npop2 = atoi( argv[arg++] ) ;
	  pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
	  for( pop =0; pop <npop2; pop++){
	    (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
	    for( i=0; i<npop2; i++){
	      if( i == pop ) arg++;
	      else {
		_mssel_argcheck( arg, argc, argv); 
		(pt->mat)[pop][i] = atof( argv[arg++] ) ;
	      }
	    }
	  }
	  for( pop = 0; pop < npop2; pop++) {
	    (pt->mat)[pop][pop] = 0.0 ;
	    for( pop2 = 0; pop2 < npop2; pop2++){
	      if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
	    }
	  }	
	}
	else {
	  _mssel_argcheck( arg, argc, argv);
	  pt->popi = atoi( argv[arg++] ) -1 ;
	  _mssel_argcheck( arg, argc, argv);
	  pt->popj = atoi( argv[arg++] ) -1 ;
	  _mssel_argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	}
	break;
      case 'j' :
	_mssel_argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	_mssel_argcheck( arg, argc, argv);
	pt->popj = atoi( argv[arg++] ) -1 ;
	break;
      default: fprintf(stderr,"e event\n");  _mssel_usage();
      }
      break;
    default: fprintf(stderr," option default\n");  _mssel_usage() ;
    }
  }
  if( (pars.mp.theta == 0.0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) ) {
    fprintf(stderr," either -s or -t or -T option must be used. \n");
    _mssel_usage();
    exit(1);
  }
  sum = 0 ;
  for( i=0; i< pars.cp.npop; i++) sum += pars.cp.config[0][i] + pars.cp.config[1][i] ;
  if( sum != pars.cp.nsam ) {
    fprintf(stderr," sum sample sizes != nsam\n");
    _mssel_usage();
    exit(1);
  }

  return 0;
}

	

/* allocates space for gametes (character strings) */
char ** _mssel_cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}

void _mssel_convertXtoADs(double *input, double *output, int n, double newmax){
  int i;
  double maxelement = -1.;
  for( i = 0; i < n; ++i)
    if(input[i] > maxelement) maxelement = input[i];

  for(i = 0; i < n; ++i)
    output[i] = input[i]/maxelement * newmax;
}

int mssel(int argc, char **argv, double ***probMatrix, int *n, int **countVector, int* densityPoints, double **ads)
{
  int i, j, k, howmany, segsites,pop, oldnfreqder , selsite, afreq, sfsMatrixCapacity = 1000, *sfsCount = NULL, *capacities = NULL ;
  char **list, **_mssel_cmatrix(), **tbsparamstrs, line[1001]  ;
  FILE *pf, *fopen() ;
  double probss, tmrca, ttot , s, age, **sfsMatrix = NULL;
  void seedit( const char * ) ;
  
  
  //void getpars( int argc, char *argv[], int *howmany )  ;
  //  static int gensam( char **list, double *probss, double *ptmrca, double *pttot, int *pselsite ) ;
  
  
  ntbs = 0 ;   /* these next few lines are for reading in parameters from a file (for each sample) */
  tbsparamstrs = (char **)malloc( argc*sizeof(char *) ) ;
  
  for( i=0; i<argc; i++) printf("%s ",argv[i]);
  for( i =0; i<argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
  for( i = 1; i<argc ; i++)
    if( strcmp( argv[i],"tbs") == 0 )  argv[i] = tbsparamstrs[ ntbs++] ;
	
  count=0;
  
  
  
  if( ntbs > 0 )  for( k=0; k<ntbs; k++)  scanf(" %s", tbsparamstrs[k] );
  int parseRes = _mssel_getpars( argc, argv, &howmany) ;   /* results are stored in global variable, pars */

  if(parseRes == 1)
    return 1;
  
  if( !pars.commandlineseedflag ) seedit( "s");
  pf = stdout ;
  
  int probMatDim = pars.cp.nsam + 1;
  
  sfsMatrix = (double**) calloc( probMatDim, sizeof(double*));
  sfsCount = (int*)calloc( probMatDim, sizeof(int));
  capacities = (int*)calloc(probMatDim, sizeof(int));
  for( i = 0; i < probMatDim; ++i){
    sfsMatrix[i] = (double*) calloc(sfsMatrixCapacity, sizeof(double));
    capacities[i] = sfsMatrixCapacity;
  }
  
  if( pars.mp.segsitesin ==  0 ) {
    list = _mssel_cmatrix(pars.cp.nsam,maxsites+1);
    posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
  }
  else {
    list = _mssel_cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
    posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
    if( pars.mp.theta > 0.0 ){
      segfac = 1.0 ;
      for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
    }
  }
  
  /*  Since the above bit is only done for count = 0, can't use tbs argument for argument of -s switch */
  
  fgets( line, 1000, pars.selfile ) ;
  sscanf(line," n: %d", &(pars.cp.nfreqder) );
  
  pars.cp.tfreqder = (double *)malloc( (unsigned)(pars.cp.nfreqder * sizeof(double) )) ;
  pars.cp.infreqder = (double **)malloc( (unsigned)(pars.cp.nfreqder * sizeof(double *))) ;
  for( i=0; i< pars.cp.nfreqder; i++)
    pars.cp.infreqder[i] = (double *)malloc( (unsigned)(pars.cp.maxpop *sizeof(double) ));
  for( i=0; i< pars.cp.nfreqder; i++){
    fscanf(pars.selfile," %lf", pars.cp.tfreqder+i );
    for(pop=0; pop < pars.cp.maxpop ; pop++){
      fscanf(pars.selfile," %lf", pars.cp.infreqder[i]+pop  );
    }
    
  }
  assert( pars.cp.nfreqder > 0 );
  /*   fprintf(stdout,"rho: %lf theta: %lf\n",pars.cp.r,pars.mp.theta); */


  
  int denPoints = 512;
  double* xpoints = calloc(denPoints+1, sizeof(double));
  double minval = xpoints[0] = 1e-8;

  /* this should be generalized for an arbitrary function
     probably it should be log or something to give more
     weight on the small distances 
  */
  /* selpos */
  double selpos = (pars.cp.selspot+1)/(double)pars.cp.nsites ;
  double maxdistance = (selpos > 0.5) ? selpos : 1. - selpos; //(double)(1.0 - 0.0)/(denPoints);
  double step = log(maxdistance/minval) /(denPoints);
  //double step = (maxdistance - minval) / denPoints ; 
  fprintf(stderr, "============== XPOINTS =================\n");
  for(i = 0; i < denPoints+1; ++i){
    xpoints[i] = exp(log(minval) + i * step);
    //xpoints[i] = minval + i * step;
    fprintf(stderr, "%.10f ", xpoints[i]);
  }

  fprintf(stderr, "=======================================\n");  
  
  while( howmany - count++ ) {
    if( count % 100 == 0)
      fprintf(stderr, "Simulation ... %d is starting... \n", count);
    
    if( (ntbs > 0) && (count >1 ) ){
      for( k=0; k<ntbs; k++){ 
	if( scanf(" %s", tbsparamstrs[k]) == EOF ){
	  if( !pars.commandlineseedflag ) seedit( "end" );
	  exit(0);
	}
      }
      _mssel_getpars( argc, argv, &howmany) ;
    }

    if( (count>1)&& ( pars.cp.nrepsfreq >1) ){
      oldnfreqder = pars.cp.nfreqder ;
      fgets(line, 1000, pars.selfile);
      fgets(line, 1000, pars.selfile);

      if( sscanf(line," n: %d", &(pars.cp.nfreqder) ) != 1 ){ 
	rewind(pars.selfile);
	do {
	  fgets( line, 1000, pars.selfile );
	}while ( line[0] == '#' );
	fgets(line, 1000, pars.selfile);
	fgets(line, 1000, pars.selfile);
	sscanf(line," n: %d", &(pars.cp.nfreqder) );
      }
      
      pars.cp.tfreqder = (double *)realloc(pars.cp.tfreqder, (unsigned)(pars.cp.nfreqder * sizeof(double) )) ;
      for( i = pars.cp.nfreqder; i < oldnfreqder;  i++) free(pars.cp.infreqder[i]);
      pars.cp.infreqder = (double **)realloc( pars.cp.infreqder,(unsigned)(pars.cp.nfreqder * sizeof(double *))) ;
      for( i= oldnfreqder; i< pars.cp.nfreqder; i++)  pars.cp.infreqder[i] = (double *)malloc( (unsigned)(pars.cp.maxpop *sizeof(double) ));

      for( i=0; i< pars.cp.nfreqder; i++){
	fscanf(pars.selfile," %lf", pars.cp.tfreqder+i );
	for(pop=0; pop < pars.cp.maxpop ; pop++)
	  fscanf(pars.selfile," %lf", pars.cp.infreqder[i]+pop  );

      }
    }
    //fprintf(pf,"\n//");
    if( ntbs >0 ){
      for(k=0; k< ntbs; k++) printf("\t%s", tbsparamstrs[k] ) ;
    }
    //printf("\n");
    segsites = _mssel_gensam( list, &probss, &tmrca, &ttot, &selsite ) ;
    if( pars.mp.timeflag ) fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
    if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
      /* if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))  */
      /*   fprintf(pf,"prob: %g\n", probss ) ; */
      
      /* if( pars.selsitepolyflag == 1 ) /\* code added by Pavlos 25/11/2015 *\/ */
      /*   fprintf(pf,"selsite: %d  %s\n", selsite, line );  /\* sites numbered 0,1,..., segsites-1 *\/ */
      
      /* fprintf(pf,"segsites: %d\n",segsites); */
      /* 	   if( segsites > 0 )	fprintf(pf,"positions: "); */
      /* 	   for( i=0; i<segsites; i++) */
      /* 		   fprintf(pf,"%6.*lf ", pars.output_precision, posit[i] ); */
      /*  fprintf(pf,"\n"); */
      if( segsites > 0 ) {
	//for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }

	
	for( i=0; i<segsites; i++){
	  for(k=afreq=0; k<pars.cp.nsam; k++) afreq += ( (list[k][i] == '1') ? 1: 0 ) ;
	  int arrayPos = afreq;
	  sfsCount[arrayPos]++;
	  if(sfsCount[arrayPos] > capacities[arrayPos]){
	    capacities[arrayPos] += capacities[arrayPos];
	    sfsMatrix[arrayPos] = (double*)realloc(sfsMatrix[arrayPos], capacities[arrayPos]*sizeof(double));
	    assert(sfsMatrix[arrayPos] != NULL);
	  }
	  sfsMatrix[arrayPos][ sfsCount[arrayPos] - 1 ] = fabs( selpos - posit[i]);
	}
      }
    }
  } // while
  
    
  for(i=0;i<probMatDim; i++) {
    qsort(sfsMatrix[i], sfsCount[i], sizeof(double), _mssel_cmpfunc);
  }

  double **probs = (double**)calloc(probMatDim, sizeof(double*));
  double **probXGivenPos = (double**)calloc(probMatDim, sizeof(double*));
  double *lads = (double*)calloc(denPoints + 1, sizeof(double)); // positions scaled, the ADs

  /* just convert the ms distances to the 0 - MAXGRID (i.e. 0 - 12 according to Nielsen) */
  _mssel_convertXtoADs(xpoints, lads, denPoints + 1, 12.0); 

  for( i = 0; i < probMatDim; ++i)
    probXGivenPos[i] = (double*)calloc(denPoints+1, sizeof(double));
    
  for(i = 0; i < probMatDim; ++i)
    probs[i] = (double*)calloc(denPoints+1, sizeof(double));
    
  double totalProb = 0.;
  double probSum = 0.;
  double *probPos = (double*)calloc( denPoints + 1 , sizeof(double));
  int totalCount = 0;
  
  for(i = 0; i < probMatDim; ++i)
    {
      kerneldensity(sfsMatrix[i], probs[i], -1, xpoints, sfsCount[i], denPoints+1, 2);
      totalCount += sfsCount[i];
    }

  /* for(i = 0; i < probMatDim; ++i) */
  /*   fprintf(pf, "sfs %d %f\n", i+1, ((double)sfsCount[i]/(double)totalCount) ); */
    
  for( i = 0; i < probMatDim; ++i)
    for(j = 0; j < denPoints+1; ++j)
      probXGivenPos[i][j] = probs[i][j] * sfsCount[i]; //(probs[i][j] * ((double)sfsCount[i]/(double)totalCount));
  

  for(j = 0; j < denPoints+1; ++j)
      for(i = 0; i < probMatDim; ++i)
	probPos[j] += probXGivenPos[i][j];
  

  for( i = 0; i < probMatDim; ++i)
    for(j = 0; j < denPoints+1; ++j)
      probXGivenPos[i][j] /= probPos[j]; //(probs[i][j] * ((double)sfsCount[i]/(double)totalCount))/probPos[j];
  

  /* for(i = 0; i < denPoints + 1; ++i) */
  /*   fprintf(pf, "%f ", xpoints[i]); */
  /* fprintf(pf, "\n"); */

  /* for(i = 0; i < probMatDim; ++i) */
  /*   { */
  /*     fprintf(pf, "%d\t", i); */
  /*     for(j = 0; j < denPoints + 1; ++j) */
  /* 	fprintf(pf, "%f ", probXGivenPos[i][j]); */
  /*     fprintf(pf, "\n"); */
  /*   } */
  
  (*probMatrix) = probXGivenPos;
  (*countVector) = sfsCount;
  (*densityPoints) = denPoints;
  (*n) = probMatDim; //pars.cp.nsam;
  (*ads) = lads;
    
  if( !pars.commandlineseedflag ) seedit( "end" );

  for(i = 0; i < argc; ++i)
    free(tbsparamstrs[i]);

  free(tbsparamstrs);
  
  free(xpoints);

  for(i = 0; i < pars.cp.nsam; ++i)
    {
      free(list[i]);
    }
  free(list);

  for(i = 0; i < probMatDim; ++i)
    free(probs[i]);
  free(probs);

  for(i = 0; i < probMatDim; ++i)
    free(sfsMatrix[i]);
  free(sfsMatrix);

  free(probPos);
  free(capacities);
}

int _mssel_pick2(n,i,j)
	int n, *i, *j;
{
  double _mssel_ran1();
  
  *i = n * _mssel_ran1() ;
  while( ( *j = n * _mssel_ran1() ) == *i );
  return(0) ;
}





/**** ordran.c  ***/






