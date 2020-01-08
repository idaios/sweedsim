/*  Link in this file for random number generation using drand48() */
/*  Thanks to Alan Rogers for suggestion of using pid.   17 Nov 2018 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

         double
ran1()
{
        double drand48();
        return( drand48() );
}               


	void seedit( char *flag )
{
	FILE *fopen(), *pfseed;
	unsigned short seedv[3], seedv2[3],  *seed48(), *pseed, tempseed ;
	int i;


	if( flag[0] == 's' ){
        time_t      currtime = time(NULL);
        unsigned long pid = (unsigned long) getpid();
        tempseed = (unsigned short)currtime^pid;
	  seedv[0] =  tempseed ;
            seedv[1] = 27011; seedv[2] = 59243; 
          seed48( seedv );   

       printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );    
	}
}


	int
commandlineseed( char **seeds)
{
	unsigned short seedv[3], *seed48();
	int i ;

	seedv[0] = atoi( seeds[0] );
	seedv[1] = atoi( seeds[1] );
	seedv[2] = atoi( seeds[2] );
	printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );    

	seed48(seedv);
	return(3);
}

