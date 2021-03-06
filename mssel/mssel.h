struct _mssel_devent {
	double time;
	int popi;
	int popj;
	double paramv;
	double **mat ;
	char detype ;
	struct _mssel_devent *nextde;
	} ;
struct _mssel_c_params {
	int npop;
	int nsam;
        int **config ;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
        int selspot ;
		int nrepsfreq;
        int nfreqder ;
        int maxpop ;
        double **infreqder ;
        double *tfreqder ;
	struct _mssel_devent *deventlist ;
	} ;
struct _mssel_m_params {
	 double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
	 } ;
struct _mssel_params { 
	struct _mssel_c_params cp;
	struct _mssel_m_params mp;
        FILE *selfile ;
        int commandlineseedflag ;
		int selsitepolyflag ;
		int output_precision;
	};
	
struct _mssel_node{
	int abv;
	int ndes;
	float time;
};

void _mssel_ordran(int n, double pbuf[]);
void _mssel_ranvec(int n, double pbuf[]);
void _mssel_order(int n, double pbuf[]);

void _mssel_biggerlist(int nsam,  char **list );
int _mssel_poisso(double u);
void _mssel_locate(int n,double beg, double len,double *ptr);
void _mssel_mnmial(int n, int nclass, double p[], int rv[]);
void _mssel_usage();
int _mssel_tdesn(struct _mssel_node *ptree, int tip, int node );
int _mssel_pick2(int n, int *i, int *j);
int _mssel_xover(int nsam,int ic, int is);
int _mssel_links(int c);

int mssel(int argc, char **argv, double ***probMatrix, int *n, int **countVector, int* densityPoints, double **ads);
