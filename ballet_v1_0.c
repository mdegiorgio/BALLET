/**********************************
 ********** BALLET v1.0 ***********
 **********************************

This software requires that you have GNU Scientific Library (GSL) installed
This software has only been tested on Linux

For speed of computation the program cycles over a predefined set of equilbrium values and recombination rates
The probality matrix for probabilities and substitutions is pre-computed prior to the start of any analysis
The probability matrix for frequency spectra is pre-computed using simulations

*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

//Variable stores probabilities of substitution and polymorphism over a grid of recombination values
//prob_poly_sub[n][0][x][r] holds probability of subsitution in sample size n that is recomination distance r at equilibrium frequency x
//prob_poly_sub[n][1][x][r] holds probability of polymorphism in sample size n that is recomination distance r at equilibrium frequency x
//Variable is global as we will use this throughout whole program
double ****prob_poly_sub = NULL;
//Variable stores probabilities of substitution and polymorphism over a grid of recombination values
//prob_spectra[n][0][x][r] holds probability of subsitution in sample size n that is recomination distance r at equilibrium frequency x
//prob_spectra[n][k][x][r] holds probability of polymorphism with ancestral frequency k in sample size n that is recomination distance r at equilibrium frequency x
//Variable is global as we will use this throughout whole program
double ****prob_spectra = NULL;
double *rate_grid = NULL; //holds vector of grid of scaled recombination rates
double *equilibria_grid = NULL; //holds vector of grid of equilibria values
int num_sites = 0;
int num_equilibria = 0;
double e_min = 0.05; // Minimum equilibrium value
double e_max = 0.951; // Maximimum equilibrium value
double e_inc = 0.05; // Equilibrium increment
int num_rates = 0;
int n_max = 0;
int n_min = 0;

double *log_fact = NULL;
double **log_binom = NULL;


struct datatype {
	double loc;
	int x;
	int n;
} *data;

struct datatype_rec {
	double loc, rate;
} *data_rec;

int num_sites_rec;
double recratemin, recratemax;

struct datatype *data=NULL;

#define NCOLUMN 3
#define LOC 0
#define X 1
#define N 2

char colname[NCOLUMN][10]={"position", "x", "n"};

#define NCOLUMN_REC 2
#define LOC_REC 0
#define RATE_REC 1

//all recombination map column headers are required.
char colname_rec[NCOLUMN_REC][10]={"position", "rate"};

void InitFact(void); 
void InitBinom(void);
void InitFact_sample_size(int); 
void InitBinom_sample_size(int);
double h(int, int, double, double, double);
double coal1(int, int, double, double, double);
double coal2(int, int, double, double, double);
double mig1(int, int, double, double, double);
double mig2(int, int, double, double, double);
gsl_vector* getDiagonal(int);
gsl_vector* getSubDiagonal(int, double, double, double);
gsl_vector* getSuperDiagonal(int, double, double, double);
double lnPower(double, double);
gsl_vector* buildConstHeightVector(int, double, double, double, gsl_vector*);
double expectedTreeHeightBalancing(int, double, double, double, double);
gsl_vector* buildConstLengthVector(int, double, double, double, gsl_vector*);
double expectedTreeLengthBalancing(int, double, double, double, double);
double ProbPoly(int, double, double, double, double, double);
void InitPolySubMatrix(double, double, double);
void ModifySpectraMatrix(char*);
void InitSpectraMatrix(double, double, double, char*);
void get_min_max_sample_sizes(void);
void get_min_max_rho(void);
void readsnps(char*);
void readrecmap(char*);
double read_divergence_time(char*);
double* read_polymorphic_sites(char*);
double* read_substitution_sites(char*);
double** read_empirical_spectra(char*, double*, double*);
double calc_ln_like_balsel(int, int, double, int, double, double);
double calc_ln_like_balsel_spectrum(int, int, double, int, double, double);
double calc_ln_like_background(int, int, double*, double*);
double calc_ln_like_background_spectrum(int, int, double**);
void scan_with_T1(char*, int, double, double, double, char*);
void scan_with_T2(char*, int, double, double, double, char*, char*);
void ComputeSpectrum(char*, int, double, double, double, int);
double estimateDivergence(double, double);
void printDivergence(char*, double);
double* calcPolyFreq(void);
double* calcSubsFreq(void);
void printFractionPolySubs(char*, double*, double*);
double** obtain_empirical_spectra(void);
void write_empirical_spectra(char*, double**);
void free_empirical_spectra(double**);
void usage(void);

int main(int argc, char *argv[])
{
	srand(time(NULL));

	int gridsize;
	int sample_size;
	int num_reps;
	double totNumSites = 0.0;
	double equilibrium_frequency;
	double rho, theta1 = 0.05, theta2 = 0.05;
	char snpfn[1000], outfn[1000], recmapfn[1000], divfn[1000], polysubfn[1000], spectfn[1000], path[10000];
	double divergence;
	double seqDiff;
	double theta;
	double *polyFreq;
	double *subsFreq;
	double **empiricalSpectra;	

	if(argc < 3) {
		usage();
	}
	
	if(strcmp(argv[1], "-inter_coal_time") == 0) {
		if(argc!=6) {
			usage();
		}
		else {
			printf("You have chosen to get mean interspecies coalescence time\n");
		}

		sprintf(snpfn, "%s", argv[2]);
		totNumSites = atof(argv[3]);
		theta = atof(argv[4]);
		sprintf(outfn, "%s", argv[5]);
		readsnps(snpfn);

		divergence = estimateDivergence(totNumSites, theta);
		printDivergence(outfn, divergence);
	}
	else if(strcmp(argv[1], "-poly_sub") == 0) {
		if(argc!=4) {
			usage();
		}
		else {
			printf("You have chosen to get proportion of polymorphic and substituted sites\n");
		}

		sprintf(snpfn, "%s", argv[2]);
		sprintf(outfn, "%s", argv[3]);
		readsnps(snpfn);

		polyFreq = calcPolyFreq();
		subsFreq = calcSubsFreq();

		printFractionPolySubs(outfn, polyFreq, subsFreq);

		free(polyFreq);
		free(subsFreq);
	}
	else if(strcmp(argv[1], "-spect") == 0) {
		if(argc!=4) {
			usage();
		}
		else {
			printf("You have chosen to get empirical spectra\n");
		}

		sprintf(snpfn, "%s", argv[2]);
		sprintf(outfn, "%s", argv[3]);
		readsnps(snpfn);

		empiricalSpectra = obtain_empirical_spectra();
		write_empirical_spectra(outfn, empiricalSpectra);
		free_empirical_spectra(empiricalSpectra);
	}
	else if(strcmp(argv[1], "-T1") == 0) {
		if(argc!=8) {
			usage();
		}
		else {
			printf("You have chosen to perform a scan using T1\n");
		}

		gridsize = atoi(argv[2]);

		if(gridsize <=0) {
			printf("windowsize should be > 0!\n");
			usage();
		}

		sprintf(divfn, "%s", argv[3]);
		sprintf(polysubfn, "%s", argv[4]);
		sprintf(snpfn, "%s", argv[5]);
		sprintf(recmapfn, "%s", argv[6]);
		sprintf(outfn, "%s", argv[7]);
		
		divergence = read_divergence_time(divfn);		
		readsnps(snpfn);
		readrecmap(recmapfn);

		InitFact();
		InitBinom();
		InitPolySubMatrix(divergence, theta1, theta2);
	
		scan_with_T1(outfn, gridsize, divergence, theta1, theta2, polysubfn);
	}
	else if(strcmp(argv[1], "-T2") == 0) {
		if(argc!=10) {
			usage();
		}
		else {
			printf("You have chosen to perform a scan using T2\n");
		}

		gridsize = atoi(argv[2]);

		if(gridsize <=0) {
			printf("windowsize should be > 0!\n");
			usage();
		}

		sprintf(divfn, "%s", argv[3]);
		sprintf(polysubfn, "%s", argv[4]);
		sprintf(spectfn, "%s", argv[5]);		
		sprintf(snpfn, "%s", argv[6]);
		sprintf(recmapfn, "%s", argv[7]);
		sprintf(path, "%s", argv[8]);
		sprintf(outfn, "%s", argv[9]);		

		divergence = read_divergence_time(divfn);		
		readsnps(snpfn);
		readrecmap(recmapfn);

		InitFact();
		InitBinom();
		InitSpectraMatrix(divergence, theta1, theta2, path);
	
		scan_with_T2(outfn, gridsize, divergence, theta1, theta2, polysubfn, spectfn);
	}
	else if(strcmp(argv[1], "-SimSpect") == 0) {
		if(argc!=8) {
			usage();
		}
		else {
			printf("You have chosen to simulate a frequency spectrum\n");
		}

		sample_size = atoi(argv[2]);
		equilibrium_frequency = atof(argv[3]);
		theta1 = atof(argv[4]);
		theta2 = atof(argv[5]);
		num_reps = atoi(argv[6]);
		sprintf(outfn, "%s", argv[7]);
		
		InitFact_sample_size(sample_size);
		InitBinom_sample_size(sample_size);
		ComputeSpectrum(outfn, sample_size, equilibrium_frequency, theta1, theta2, num_reps);
	}
	else {
		usage();
	}

	return 0;
}




void InitFact(void) 
{	
	int n = 0;
	
	log_fact = (double*)malloc((n_max+1)*sizeof(double));
	log_fact[0] = 0.0;
	for(n = 1; n <= n_max; n++) {
		log_fact[n] = log_fact[n-1] + log((double)n);
	}
}

void InitBinom(void)
{
	int n = 0;
	int k = 0;

	printf("Initializing binomial coefficients\n");

	log_binom = (double**)malloc((n_max+1)*sizeof(double*));
	for(n = n_min; n <= n_max; n++) {
		log_binom[n] = (double*)malloc((n_max+1)*sizeof(double));
		
		for(k = 0; k <= n; k++) {
			log_binom[n][k] = log_fact[n] - log_fact[n-k] - log_fact[k];
		}
	}
}

void InitFact_sample_size(int sample_size) 
{	
	int n = 0;
	
	log_fact = (double*)malloc((sample_size+1)*sizeof(double));
	log_fact[0] = 0.0;
	for(n = 1; n <= sample_size; n++) {
		log_fact[n] = log_fact[n-1] + log((double)n);
	}
}

void InitBinom_sample_size(int sample_size)
{
	int n = 0;
	int k = 0;

	printf("Initializing binomial coefficients\n");

	log_binom = (double**)malloc((sample_size+1)*sizeof(double*));
	for(n = 1; n <= sample_size; n++) {
		log_binom[n] = (double*)malloc((sample_size+1)*sizeof(double));
		
		for(k = 0; k <= n; k++) {
			log_binom[n][k] = log_fact[n] - log_fact[n-k] - log_fact[k];
		}
	}
}

double h(int i, int j, double x, double b1, double b2)
{
	if((i <= 1) && (j <= 1)) {
		return ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
	else if(i <= 1) {
		return ((j*(j-1)/2)/(1-x)) + ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
	else if(j <= 1) {
		return ((i*(i-1)/2)/x) + ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
	else {
		return ((i*(i-1)/2)/x) + ((j*(j-1)/2)/(1-x)) + ((j*b1*x)/(1-x)) + ((i*b2*(1-x))/x);
	}
}

double coal1(int i, int j, double x, double b1, double b2)
{
	if(i <= 1) {
		return 0.0;
	}
	else {
		return (i*(i-1)/2)/(x*h(i,j,x,b1,b2));
	}
}

double coal2(int i, int j, double x, double b1, double b2)
{
	if(j <= 1) {
		return 0.0;
	}
	else {
		return (j*(j-1)/2)/((1-x)*h(i,j,x,b1,b2));
	}
}

double mig1(int i, int j, double x, double b1, double b2)
{
	return (i*b2*(1-x)) / (x * h(i, j, x, b1, b2));
}

double mig2(int i, int j, double x, double b1, double b2)
{
	return (j*b1*x) / ((1-x)*h(i, j, x, b1, b2));
}

gsl_vector* getDiagonal(int n)
{
	int i = 0;
	gsl_vector *diag = gsl_vector_alloc(n+1);

	for(i = 0; i < n+1; i++) {
		gsl_vector_set(diag, i, 1.0);
	}

	return diag;
}

gsl_vector* getSubDiagonal(int n, double x, double b1, double b2)
{
	int i = 0;
	gsl_vector *subDiag = gsl_vector_alloc(n);

	for(i = 0; i < n; i++) {
		gsl_vector_set(subDiag, i, -mig1(i+1, n - (i+1), x, b1, b2));
	}

	return subDiag;
}

gsl_vector* getSuperDiagonal(int n, double x, double b1, double b2)
{
	int i = 0;
	gsl_vector *superDiag = gsl_vector_alloc(n);

	for(i = 0; i < n; i++) {
		gsl_vector_set(superDiag, i, -mig2(i, n - i, x, b1, b2));
	}

	return superDiag;
}

double lnPower(double x, double n)
{
	if(x <= 0.0) {
		printf("ERROR: x <= 0.0 in lnPower()\n");//// MAKE THIS MESSAGE A LITTLE BETTER
		exit(-1);
	}

	return ((double)n)*log((double)x);
}

gsl_vector* buildConstHeightVector(int n, double x, double b1, double b2, gsl_vector* L)
{
	int k;
	gsl_vector *result = gsl_vector_alloc(n+1);

	// COMPUTES result[0] = (1 / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * L[0]
	gsl_vector_set(result, 0, (1 / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * gsl_vector_get(L,0));
	// COMPUTES result[n] = (1 / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * L[n-1]
	gsl_vector_set(result, n, (1 / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * gsl_vector_get(L,n-1));

	for(k = 1; k <= n -1; k++) {
		// COMPUTES result[k] = (1 / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * L[k-1] + coal2(k, n - k, x, b1, b2) * L[k];
		gsl_vector_set(result, k, (1 / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * gsl_vector_get(L,k-1) + coal2(k, n - k, x, b1, b2) * gsl_vector_get(L,k));
	}

	return result;
}

double expectedTreeHeightBalancing(int n, double x, double theta1, double theta2, double R)
{
	int k;
	int val = 0;
	double treeHeight = 0;
	double b1 = theta1 + R*(1 - x);
	double b2 = theta2 + R*x;
	gsl_vector *diag;
	gsl_vector *subDiag;
	gsl_vector *superDiag;
	gsl_vector *u;
	gsl_vector *v;

	v = gsl_vector_alloc(2);
	gsl_vector_set(v,0,0);
	gsl_vector_set(v,1,0);

	for(k = 2; k <= n; k++) {
		u = buildConstHeightVector(k, x, b1, b2, v);
		gsl_vector_free(v);
		v = gsl_vector_alloc(k+1);

		diag = getDiagonal(k);
		subDiag = getSubDiagonal(k, x, b1, b2);
		superDiag = getSuperDiagonal(k, x, b1, b2);

		val = gsl_linalg_solve_tridiag(diag, superDiag, subDiag, u, v);

		gsl_vector_free(u);
		gsl_vector_free(diag);
		gsl_vector_free(subDiag);
		gsl_vector_free(superDiag);
	}

	for(k = 0; k <= n; k++) {
		treeHeight = treeHeight + gsl_vector_get(v,k) * exp(log_binom[n][k] + lnPower(x, k) + lnPower(1-x, n - k));
	}

	gsl_vector_free(v);

	return treeHeight;
}

gsl_vector* buildConstLengthVector(int n, double x, double b1, double b2, gsl_vector* L)
{
	int k;
	gsl_vector *result = gsl_vector_alloc(n+1);

	// COMPUTES result[0] = (n / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * L[0]
	gsl_vector_set(result, 0, (n / h(0, n, x, b1, b2))  + coal2(0, n, x, b1, b2) * gsl_vector_get(L,0));
	// COMPUTES result[n] = (n / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * L[n-1]
	gsl_vector_set(result, n, (n / h(n, 0, x, b1, b2))  + coal1(n, 0, x, b1, b2) * gsl_vector_get(L,n-1));

	for(k = 1; k <= n -1; k++) {
		// COMPUTES result[k] = (n / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * L[k-1] + coal2(k, n - k, x, b1, b2) * L[k];
		gsl_vector_set(result, k, (n / h(k, n - k, x, b1, b2)) + coal1(k, n - k, x, b1, b2) * gsl_vector_get(L,k-1) + coal2(k, n - k, x, b1, b2) * gsl_vector_get(L,k));
	}

	return result;
}

double expectedTreeLengthBalancing(int n, double x, double theta1, double theta2, double R)
{
	int k;
	int val = 0;
	double treeLen = 0;
	double b1 = theta1 + R*(1 - x);
	double b2 = theta2 + R*x;
	gsl_vector *diag;
	gsl_vector *subDiag;
	gsl_vector *superDiag;
	gsl_vector *u;
	gsl_vector *v;

	v = gsl_vector_alloc(2);
	gsl_vector_set(v,0,0);
	gsl_vector_set(v,1,0);

	for(k = 2; k <= n; k++) {
		u = buildConstLengthVector(k, x, b1, b2, v);
		gsl_vector_free(v);
		v = gsl_vector_alloc(k+1);

		diag = getDiagonal(k);
		subDiag = getSubDiagonal(k, x, b1, b2);
		superDiag = getSuperDiagonal(k, x, b1, b2);

		val = gsl_linalg_solve_tridiag(diag, superDiag, subDiag, u, v);

		gsl_vector_free(u);
		gsl_vector_free(diag);
		gsl_vector_free(subDiag);
		gsl_vector_free(superDiag);
	}

	for(k = 0; k <= n; k++) {
		treeLen = treeLen + gsl_vector_get(v,k) * exp(log_binom[n][k] + lnPower(x, k) + lnPower(1-x, n - k));
	}

	gsl_vector_free(v);

	return treeLen;
}

double ProbPoly(int n, double divergence, double x, double theta1, double theta2, double R)
{
	double treeHeight = 0.0;
	double treeLen = 0.0;

	treeLen = expectedTreeLengthBalancing(n, x, theta1, theta2, R);
	treeHeight = expectedTreeHeightBalancing(n, x, theta1, theta2, R);

	if(divergence <= treeHeight) {
		printf("ERROR: Expected tree height (%lf) > ingroup-outgroup expected coalescence time (%lf) in selection likelihood function\n", treeHeight, divergence);
		printf("You need to choose a larger theta values (theta1,theta2)=(%lf,%lf)\n", theta1, theta2);
		exit(-1);
	}

	return exp(log(treeLen) - log(treeLen + 2.0*divergence - treeHeight));
}



void InitPolySubMatrix(double divergence, double theta1, double theta2) 
{	
	int n = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int count = 0;
	double x = 0.0;	
	double r = 0.0;	

	printf("Initializing matrix of polymorphism and substitution probabilities\n");

	// Identify the number of equilibrium values
	for(x = e_min; x <= e_max; x = x + e_inc) {
		num_equilibria++;
	}

	// Initialize equilibria grid
	equilibria_grid = (double*)malloc(num_equilibria*sizeof(double));
	for(x = e_min; x <= e_max; x = x + e_inc) {
		equilibria_grid[count] = x;	
		count++;
	}

	// Identify the number of recombination values
	for(r = 0.0; r <= 1.0; r = r + 0.001) {
		num_rates++;
	}
	for(r = 1.1; r <= 10.0; r = r + 0.1) {
		num_rates++;
	}
	for(r = 20.0; r <= 100.0; r = r + 10.0) {
		num_rates++;
	}
		
	// Initialize recombination rate grid
	rate_grid = (double*)malloc(num_rates*sizeof(double));
	count = 0;	
	for(r = 0.0; r <= 1.0; r = r + 0.001) {
		rate_grid[count] = r;
		count++;
	}
	for(r = 1.1; r <= 10.0; r = r + 0.1) {
		rate_grid[count] = r;
		count++;
	}
	for(r = 20.0; r <= 100.0; r = r + 10.0) {
		rate_grid[count] = r;
		count++;
	}

	// Initialize probability matrix
	prob_poly_sub = (double****)malloc((n_max+1)*sizeof(double***));
	for(n = n_min; n <= n_max; n++) {
		prob_poly_sub[n] = (double***)malloc(2*sizeof(double**));

		for(j = 0; j < 2; j++) {
			prob_poly_sub[n][j] = (double**)malloc(num_rates*sizeof(double*));
			
			for(k = 0; k < num_equilibria; k++) {
				prob_poly_sub[n][j][k] = (double*)malloc(num_rates*sizeof(double));

				for(l = 0; l < num_rates; l++) {
					prob_poly_sub[n][j][k][l] = 0.0;   //////////// REPLACE THIS WITH AN ACTUAL PROBABILITY UNDER THE MODEL
				}
			}
		}	

		for(k = 0; k < num_equilibria; k++) {
			for(l = 0; l < num_rates; l++) {
				prob_poly_sub[n][1][k][l] = ProbPoly(n, divergence, equilibria_grid[k], theta1, theta2, rate_grid[l]);
				prob_poly_sub[n][0][k][l] = 1.0 - prob_poly_sub[n][1][k][l];
			}
		}
	}
}

void ModifySpectraMatrix(char *path)
{
	FILE *infile;
	int n, k, l, j;	
	double total;
	char filename[1000];
	double *spectrum = NULL;
	
	spectrum = (double*)malloc(n_max * sizeof(double));
	
	// *******SHOULD MAYBE CHANGE TO INCREMEMENTS OF 1 RATHER THAN 2
	for(n = n_min; n <= n_max; n = n + 2) { 
		for(k = 0; k < num_equilibria; k++) {
			sprintf(filename, "%s/n%d/spectrum_n%d_x%d.txt", path, n, n, 5*k + 5);
		//	printf("Opening %s\n", filename);
						
			infile = fopen(filename, "r");			
	
			for(l = 0; l < num_rates; l++) {
				for(j = 0; j < n; j++) {
					fscanf(infile, "%lf", &spectrum[j]); // Read a line from the file
				}

				total = 0.0;
				for(j = 1; j < n; j++) {
					total = total + spectrum[j];	
				}
				
				for(j = 1; j < n; j++) {
					spectrum[j] = spectrum[j] / total;	
				}
				
				// Convert from derived to ancestral and set to spectum variable
				for(j = 1; j < n; j++) {
					prob_spectra[n][j][k][l] = (1.0 - prob_spectra[n][0][k][l]) * spectrum[n-j];   
				}				
			}
		}	
	}

	free(spectrum);
}


void InitSpectraMatrix(double divergence, double theta1, double theta2, char *path) 
{	
	int n = 0;
	int j = 0;
	int k = 0;
	int l = 0;
	int count = 0;
	double x = 0.0;	
	double r = 0.0;	

	printf("Initializing matrix of spectra probabilities\n");

	// Identify the number of equilibrium values
	for(x = e_min; x <= e_max; x = x + e_inc) {
		num_equilibria++;
	}

	// Initialize equilibria grid
	equilibria_grid = (double*)malloc(num_equilibria*sizeof(double));
	for(x = e_min; x <= e_max; x = x + e_inc) {
		equilibria_grid[count] = x;	
		count++;
	}

	// Identify the number of recombination values
	for(r = 0.0; r <= 1.0; r = r + 0.001) {
		num_rates++;
	}
	for(r = 1.1; r <= 10.0; r = r + 0.1) {
		num_rates++;
	}
	for(r = 20.0; r <= 100.0; r = r + 10.0) {
		num_rates++;
	}
		
	// Initialize recombination rate grid
	rate_grid = (double*)malloc(num_rates*sizeof(double));
	count = 0;	
	for(r = 0.0; r <= 1.0; r = r + 0.001) {
		rate_grid[count] = r;
		count++;
	}
	for(r = 1.1; r <= 10.0; r = r + 0.1) {
		rate_grid[count] = r;
		count++;
	}
	for(r = 20.0; r <= 100.0; r = r + 10.0) {
		rate_grid[count] = r;
		count++;
	}

	// Initialize spectra matrix
	prob_spectra = (double****)malloc((n_max+1)*sizeof(double***));
	for(n = n_min; n <= n_max; n++) {
		prob_spectra[n] = (double***)malloc(n*sizeof(double**));

		for(j = 0; j < n; j++) {
			prob_spectra[n][j] = (double**)malloc(num_rates*sizeof(double*));
			
			for(k = 0; k < num_equilibria; k++) {
				prob_spectra[n][j][k] = (double*)malloc(num_rates*sizeof(double));

				for(l = 0; l < num_rates; l++) {
					prob_spectra[n][j][k][l] = 0.0;   //////////// REPLACE THIS WITH AN ACTUAL PROBABILITY UNDER THE MODEL
				}
			}
		}	

		for(k = 0; k < num_equilibria; k++) {
			for(l = 0; l < num_rates; l++) {
				prob_spectra[n][0][k][l] = 1.0 - ProbPoly(n, divergence, equilibria_grid[k], theta1, theta2, rate_grid[l]);
			}
		}
	}
	ModifySpectraMatrix(path);
}

void get_min_max_sample_sizes(void)
{
	int i;

	if(num_sites==0) {
		n_min = 0;
		n_max = 0;
		return;
	}

	n_min = data[0].n;
	n_max = data[0].n;

	for(i = 0; i < num_sites; i++) {
		if(data[i].n > n_max) {
			n_max = data[i].n;
		}

		if(data[i].n < n_min) {
			n_min = data[i].n;
		}
	}
}


void get_min_max_rho(void)
{
	int i;

	if(num_sites_rec==0) {
		recratemin = 0.0;
		recratemax = 0.0;
		return;
	}

	recratemin = data_rec[0].rate;
	recratemax = data_rec[0].rate;

	for(i = 0; i < num_sites_rec; i++) {
		if(data_rec[i].rate > recratemax) {
			recratemax = data_rec[i].rate;
		}

		if(data_rec[i].rate < recratemin) {
			recratemin = data_rec[i].rate;
		}
	}
}

void readsnps(char *infn)
{
	FILE *infile=fopen(infn, "r");
	int colpos[NCOLUMN], pos=0, col=0, i, j;
	char c, str[1000];

	printf("Reading genetic variation file %s\n", infn);

	if(num_sites > 0) {
		free(data);
	}

	num_sites=0;
	c=fgetc(infile);

	for(i=0; i<NCOLUMN; i++) {
		colpos[i]=-1;
	}

	while(c!='\n' && c!=EOF) {
		str[0]=c;
		pos=1;

		while('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF) {
			str[pos++]=c;
		}

		str[pos]='\0';

		for(i=0; i<NCOLUMN; i++) {
			if(strcmp(str, colname[i])==0) {
				colpos[i]=col;
				break;
			}
		}

		col++;

		if(c=='\n' || c==EOF) {
			break;
		}

		c=fgetc(infile);
	}

	if(colpos[LOC]<0 || colpos[X]<0 || colpos[N]<0) {
		fprintf(stderr, "readsnps: infile should have columns named position, x, and n\n");
		exit(-1);
	}

	while(EOF!=(c=fgetc(infile))) {
		if(c=='\n') {
			num_sites++;
		}
	}

	fclose(infile);
	infile=fopen(infn, "r");

	while('\n'!=(c=fgetc(infile)) && c!=EOF) {
		;
	}

	data = malloc(num_sites*sizeof(struct datatype));

	for(i=0; i<num_sites; i++) {
		for(j=0; j<col; j++) {
			if(colpos[LOC]==j) {
				if(EOF==fscanf(infile, "%lf", &data[i].loc)) {
					printf("Could not read %s\n", infn);
					exit(-1);
				}
			}
			else if(colpos[X]==j) {
				if(EOF==fscanf(infile, "%i", &data[i].x)) {
					printf("Could not read %s\n", infn);
					exit(-1);
				}
			}
			else if(colpos[N]==j) {
				if(EOF==fscanf(infile, "%i", &data[i].n)) {
					printf("Could not read %s\n", infn);
					exit(-1);
				}
			}
			else {
				while('\t'!=(c=fgetc(infile)) && c!='\n' && c!=EOF) {
					;
				}

				if(c=='\n' || c==EOF) {
					if(j!=col-1) {
						printf("Could not read %s\n", infn);
						exit(-1);
					}
				}
			}
		}
	}

	fclose(infile);

	get_min_max_sample_sizes();

	printf("\tDone reading SNPs (num_sites, min_sample_size, max_sample_size) = (%d, %d, %d)\n", num_sites, n_min, n_max);
}

void readrecmap(char *infn)
{
	FILE *infile = fopen(infn, "r");
	int colpos[NCOLUMN_REC], pos=0, col=0, i, j;
	char c, str[1000];

	printf("Reading recombination file %s\n", infn);

	if(num_sites_rec > 0) {
		free(data_rec);
	}

	num_sites_rec = 0;
	c = fgetc(infile);

	for(i = 0; i < NCOLUMN_REC; i++) {
		colpos[i] = -1;
	}

	while(c!='\n' && c!=EOF) {
		str[0] = c;
		pos = 1;

		while('\t' != (c = fgetc(infile)) && c != '\n' && c != EOF) {
			str[pos++] = c;
		}

		str[pos] = '\0';

		for(i = 0; i < NCOLUMN_REC; i++) {
			if(strcmp(str, colname_rec[i]) == 0) {
				colpos[i] = col;
				break;
			}
		}

		col++;

		if(c == '\n' || c == EOF) {
			break;
		}

		c = fgetc(infile);
	}

	if(colpos[LOC_REC] < 0 || colpos[RATE_REC] < 0) {
		fprintf(stderr, "readrecmap: infile should have columns named position and rate\n");
		exit(-1);
	}

	while(EOF != (c = fgetc(infile))) {
		if(c == '\n') {
			num_sites_rec++;
		}
	}

	fclose(infile);
	infile = fopen(infn, "r");

	while('\n' != (c = fgetc(infile)) && c != EOF) {
		;
	}

	data_rec = malloc(num_sites_rec * sizeof(struct datatype_rec));

	for(i = 0; i < num_sites; i++) {
		for(j = 0; j < col; j++) {
			if(colpos[LOC_REC] == j) {
				if(EOF == fscanf(infile, "%lf", &data_rec[i].loc)) {
					printf("Could not read %s\n", infn);
					exit(-1);
				}
			}
			else if(colpos[RATE_REC] == j) {
				if(EOF == fscanf(infile, "%le", &data_rec[i].rate)) {
					printf("Could not read %s\n", infn);
					exit(-1);
				}

				// Constructing cumulative rate
				data_rec[i].rate = data_rec[i].rate + data_rec[i - 1].rate;
			}
			else {
				while('\t' != (c = fgetc(infile)) && c != '\n' && c != EOF) {
					;
				}

				if(c == '\n' || c == EOF) {
					if(j != col - 1) {
						printf("Could not read %s\n", infn);
						exit(-1);
					}
				}
			}
		}
	}

	fclose(infile);

	get_min_max_rho();

	printf("\tDone reading recombination rates cumulative_rho = %lf\n", recratemax);
}

double read_divergence_time(char *infn)
{
	FILE *infile;
	double divergence;

	infile=fopen(infn, "r");

	fscanf(infile, "%lf", &divergence);

	fclose(infile);

	return divergence;
}


double* read_polymorphic_sites(char *polysubfn)
{
	FILE *infile;
	int i;
	double *freq;
	int dummy1;
	double dummy2;

	freq = (double *)malloc((n_max+1) * sizeof(double));

	infile=fopen(polysubfn, "r");

	for(i = 0; i <= n_max; i++) {
		fscanf(infile, "%d", &dummy1);
		fscanf(infile, "%lf", &freq[i]);
		fscanf(infile, "%lf", &dummy2);
	}

	fclose(infile);

	return freq;

}

double* read_substitution_sites(char *polysubfn)
{
	FILE *infile;
	int i;
	double *freq;
	int dummy1;
	double dummy2;

	freq = (double *)malloc((n_max+1) * sizeof(double));

	infile=fopen(polysubfn, "r");

	for(i = 0; i <= n_max; i++) {
		fscanf(infile, "%d", &dummy1);
		fscanf(infile, "%lf", &dummy2);
		fscanf(infile, "%lf", &freq[i]);
	}

	fclose(infile);

	return freq;

}

double** read_empirical_spectra(char *spectfn, double *polyFreq, double *subsFreq)
{
	FILE* infile;
	int i, j, minSample, maxSample;
	double total = 0.0;
	double **spectra;

	spectra = (double **)malloc((n_max+1) * sizeof(double *));

	for(i = 0; i <= n_max; i++) {
		spectra[i] = (double *)malloc((i+1) * sizeof(double));
	}

	printf("reading %s\n", spectfn);
	infile = fopen(spectfn, "r");
	
	fscanf(infile, "%d", &minSample);
	fscanf(infile, "%d", &maxSample);
	for(i = minSample; i <= maxSample; i++) {
		total = 0.0;
		for(j = 1; j <= i - 1; j++) {
			fscanf(infile, "%lf", &spectra[i][i - j]);
			total = total + spectra[i][i - j];		
		}
		
		spectra[i][0] = subsFreq[i];
		for(j = 1; j <= i - 1; j++) {
			spectra[i][j] = polyFreq[i] * (spectra[i][j] / total);
		}
	}
	
	return spectra;
}

double calc_ln_like_balsel(int marker, int windowRadius, double divergence, int x, double theta1, double theta2)
{
	double ln_like = 0.0;
	double recombinationRate = 0.0;
	int i = 0;
	int r = 0; // VARIABLE FOR RECOMBINATION RATE INDEX
	double interval_factor = 0.0;

	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
		recombinationRate = fabs(data_rec[i].rate - data_rec[marker].rate);
		
		// If the recombination rate is higher than the number of rates in the file, just use the spectrum for the largest rate.
		// This should be fine because the largest rate should be large enough such that you are close to neutrality. 
		if(rate_grid[num_rates-1] < recombinationRate) { 
			r = num_rates - 1;
		}

		// Find recombination rate index that fits recombination rate
		for(r = 0; r < num_rates-1; r++) {
			// If found index
			if(rate_grid[r] <= recombinationRate && recombinationRate <= rate_grid[r + 1]) {
				break;				
			}
		}

		interval_factor = (recombinationRate - rate_grid[r]) / (rate_grid[r + 1] - rate_grid[r]);

		if(data[i].x == 0) {
			ln_like = ln_like + log(prob_poly_sub[data[i].n][0][x][r] + interval_factor * (prob_poly_sub[data[i].n][0][x][r + 1] - prob_poly_sub[data[i].n][0][x][r]));
			
		}
		else if(data[i].x > 0 && data[i].x < data[i].n) {
			ln_like = ln_like + log(prob_poly_sub[data[i].n][1][x][r] + interval_factor * (prob_poly_sub[data[i].n][1][x][r + 1] - prob_poly_sub[data[i].n][1][x][r]));
		}
	}

	return ln_like;
}

double calc_ln_like_balsel_spectrum(int marker, int windowRadius, double divergence, int x, double theta1, double theta2)
{
	double ln_like = 0.0;
	double recombinationRate = 0.0;
	int i = 0;
	int r = 0; // VARIABLE FOR RECOMBINATION RATE INDEX
	double interval_factor = 0.0;

	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
		recombinationRate = fabs(data_rec[i].rate - data_rec[marker].rate);
		
		// If the recombination rate is higher than the number of rates in the file, just use the spectrum for the largest rate.
		// This should be fine because the largest rate should be large enough such that you are close to neutrality. 
		if(rate_grid[num_rates-1] < recombinationRate) { 
			r = num_rates - 1;
		}

		// Find recombination rate index that fits recombination rate
		for(r = 0; r < num_rates-1; r++) {
			// If found index
			if(rate_grid[r] <= recombinationRate && recombinationRate <= rate_grid[r + 1]) {
				break;				
			}
		}

		interval_factor = (recombinationRate - rate_grid[r]) / (rate_grid[r + 1] - rate_grid[r]);

		if(data[i].x == 0) {
			ln_like = ln_like + log(prob_spectra[data[i].n][0][x][r] + interval_factor * (prob_spectra[data[i].n][0][x][r + 1] - prob_spectra[data[i].n][0][x][r]));
			
		}
		else if(data[i].x > 0 && data[i].x < data[i].n) {
			ln_like = ln_like + log(prob_spectra[data[i].n][data[i].x][x][r] + interval_factor * (prob_spectra[data[i].n][data[i].x][x][r + 1] - prob_spectra[data[i].n][data[i].x][x][r]));
		}
	}

	return ln_like;
}

double calc_ln_like_background(int marker, int windowRadius, double *polyFreq, double *subsFreq)
{
	double ln_like = 0.0;
	int i = 0;

	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
		if(data[i].x == 0) {
			ln_like = ln_like + log(subsFreq[data[i].n]);
		}
		else if(data[i].x > 0 && data[i].x < data[i].n) {
			ln_like = ln_like + log(polyFreq[data[i].n]);
		}
	}

	return ln_like;
}

double calc_ln_like_background_spectrum(int marker, int windowRadius, double **spectra)
{
	double ln_like = 0.0;
	int i = 0;

	for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
		if(data[i].x == 0) {
			ln_like = ln_like + log(spectra[data[i].n][0]);
		}
		else if(data[i].x > 0 && data[i].x < data[i].n) {
			ln_like = ln_like + log(spectra[data[i].n][data[i].x]);
		}
	}

	if(ln_like != ln_like) { // If NaN
		printf("Found a NaN at marker %d with spectra\n", marker);
		for(i = marker - windowRadius; i <= marker + windowRadius; i++) {
			printf("\tn=%d\tx=%d\tspectrum=%lf\tlogSpectrum=%lf\n", data[i].n, data[i].x, spectra[data[i].n][data[i].x], log(spectra[data[i].n][data[i].x]));		
		}
		exit(-1);
	}

	return ln_like;
}

void scan_with_T1(char *outfn, int windowRadius, double divergence, double theta1, double theta2, char *polysubfn)
{
	int startPos = windowRadius;
	int endPos = num_sites - windowRadius - 1;
	int i = 0;
	FILE* outfile;
	double ln_like_balsel, ln_like_background, ln_lr;
	double ln_like_temp;	
	int x_max, x;	
	double *polyFreq;
	double *subsFreq;
	
	polyFreq = read_polymorphic_sites(polysubfn);
	subsFreq = read_substitution_sites(polysubfn);

	outfile=fopen(outfn, "w");

	printf("Performing scan using T1 and writing results to %s...\n", outfn);

	for(i = startPos; i <= endPos; i++) {
		x_max = equilibria_grid[0];
		ln_like_balsel = calc_ln_like_balsel(i, windowRadius, divergence, x_max, theta1, theta2);

		for(x = 1; x < num_equilibria; x++) {
			ln_like_temp = calc_ln_like_balsel(i, windowRadius, divergence, x, theta1, theta2);		
			
			if(ln_like_balsel < ln_like_temp) {
				ln_like_balsel = ln_like_temp;
				x_max = x;
			}		
		}

		ln_like_background = calc_ln_like_background(i, windowRadius, polyFreq, subsFreq);
		ln_lr = 2*(ln_like_balsel - ln_like_background);
		fprintf(outfile, "%lf\t%lf\n", data[i].loc, ln_lr);
	}

	fclose(outfile);

	free(polyFreq);
	free(subsFreq);
}

void scan_with_T2(char *outfn, int windowRadius, double divergence, double theta1, double theta2, char *polysubfn, char* spectfn)
{
	int startPos = windowRadius;
	int endPos = num_sites - windowRadius - 1;
	int i = 0;
	FILE* outfile;
	double ln_like_balsel, ln_like_background, ln_lr;
	double ln_like_temp;	
	int x_max, x;	
	double *polyFreq;
	double *subsFreq;
	double ***spectra;
	double **empiricalSpectra;
	
	polyFreq = read_polymorphic_sites(polysubfn);
	subsFreq = read_substitution_sites(polysubfn);
	empiricalSpectra = read_empirical_spectra(spectfn, polyFreq, subsFreq);
	
	free(polyFreq);
	free(subsFreq);

	outfile=fopen(outfn, "w");

	printf("Performing scan using T2 and writing results to %s...\n", outfn);

	for(i = startPos; i <= endPos; i++) {
		x_max = equilibria_grid[0];
		ln_like_balsel = calc_ln_like_balsel_spectrum(i, windowRadius, divergence, x_max, theta1, theta2);

		for(x = 1; x < num_equilibria; x++) {
			ln_like_temp = calc_ln_like_balsel_spectrum(i, windowRadius, divergence, x, theta1, theta2);		
			
			if(ln_like_balsel < ln_like_temp) {
				ln_like_balsel = ln_like_temp;
				x_max = x;
			}		
		}
		
		ln_like_background = calc_ln_like_background_spectrum(i, windowRadius, empiricalSpectra);
		ln_lr = 2*(ln_like_balsel - ln_like_background);
		fprintf(outfile, "%lf\t%lf\n", data[i].loc, ln_lr);
	}

	fclose(outfile);
}

void ComputeSpectrum(char *outfn, int n, double x, double theta1, double theta2, int num_reps)
{
	int k = 0;
	int i = 0;
	int n1 = 0;
	int n2 = 0;
	int rep = 0;
	int choice1 = 0;
	int choice2 = 0;
	int count = 0;
	int start_count = 0;
	double R = 0.0;
	double total = 0.0;
	double rand_num = 0.0;
	double sum_choice = 0.0;
	double v_h = 0.0;
	double v_c1 = 0.0;
	double v_c2 = 0.0;
	double v_m1 = 0.0;
	double v_m2 = 0.0; 
	double beta1 = 0.0;
	double beta2 = 0.0;
	double *times = NULL;
	double *tot_times = NULL;
	double *v1 = NULL;
	double *v2 = NULL;
	FILE* outfile;
	char command[1000];
	char line[1000000];

	// Identify the number of recombination values
	for(R = 0.0; R <= 1.0; R = R + 0.001) {
		num_rates++;
	}
	for(R = 1.1; R <= 10.0; R = R + 0.1) {
		num_rates++;
	}
	for(R = 20.0; R <= 100.0; R = R + 10.0) {
		num_rates++;
	}
		
	// Initialize recombination rate grid
	rate_grid = (double*)malloc(num_rates*sizeof(double));
	count = 0;	
	for(R = 0.0; R <= 1.0; R = R + 0.001) {
		rate_grid[count] = R;
		count++;
	}
	for(R = 1.1; R <= 10.0; R = R + 0.1) {
		rate_grid[count] = R;
		count++;
	}
	for(R = 20.0; R <= 100.0; R = R + 10.0) {
		rate_grid[count] = R;
		count++;
	}

	times = (double*)malloc((n+1)*sizeof(double));
	tot_times = (double*)malloc((n+1)*sizeof(double));
	v1 = (double*)malloc((n+1)*sizeof(double));
	v2 = (double*)malloc((n+1)*sizeof(double));

	outfile=fopen(outfn, "r");
	if(outfile == NULL) { // If the simulated spectrum file doesn't exist
		outfile=fopen(outfn, "w");
	}
	else { // If the simulated spectrum file exists
		while(fgets(line, sizeof(line), outfile)) {
			start_count++;
			//printf("count=%d %s", start_count, line);
		}
		fclose(outfile);
	
		if(start_count > 0) {
			start_count--;
			sprintf(command, "cp %s temp_%s", outfn, outfn);
			system(command);//printf("command is %s\n", command);		
			sprintf(command, "head -%d temp_%s > %s", start_count, outfn, outfn);
			system(command);//printf("command is %s\n", command);	
			sprintf(command, "rm -f temp_%s", outfn);
			system(command);//printf("command is %s\n", command);	
			outfile=fopen(outfn, "a");	
		}
		else {
			outfile=fopen(outfn, "w");
		}
	}
	

	for(count = start_count; count < num_rates; count++) {
		for(i = 1; i <= n-1; i++) {
			tot_times[i] = 0.0;			
		}

		R = rate_grid[count];
		beta1 = theta1 + R*(1.0-x);
		beta2 = theta2 + R*x;

		printf("\tGetting spectrum for (n,x,rho,beta1,beta2) = (%d,%lf,%lf,%lf,%lf) with rep %d/%d ", n, x, R, beta1, beta2, count+1, num_rates);

		for(rep = 1; rep <= num_reps; rep++) {
			// Try possible k values
			for(k = 0; k <= n; k++) {
				// Initialize vector of states
				for(i = 1; i <= n; i++) {
					times[i] = 0.0;
					v1[i] = 0.0;
					v2[i] = 0.0;
				}
				v1[1] = k;
				v2[1] = n - k;

				n1 = k;
				n2 = n- k;

				while(1) {
					v_h = n1*(n1-1.0)/(2.0*x) + n2*(n2-1.0)/(2.0*(1.0-x)) + n2*beta1*(x/(1.0-x)) + n1*beta2*((1.0-x)/x);
					v_c1 = n1*(n1-1.0)/(2.0*x*v_h);
					v_c2 = n2*(n2-1.0)/(2.0*(1.0-x)*v_h);
					v_m1 = n1*beta2*((1.0-x)/x)/v_h;
					v_m2 = n2*beta1*(x/(1.0-x))/v_h;
				
					rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
					while(rand_num == 0.0){ // Guard against log(0)
						printf("Guarding against log(0)\n");
						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
					}
					rand_num = -(1.0/v_h)*log(rand_num); // Genereate exponential RV

					// Compute times
					for(i = 1; i <= n-1; i++) {
						times[i] = times[i] + ((double)v1[i]+v2[i]) * rand_num;
					}

					rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
				
					if(rand_num <= v_c1) {
						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
						choice1 = -1;
						choice2 = -1;
						sum_choice = 0.0;
						// Look or non-empty cells
						for(i = 1; i <= n; i++) {
							sum_choice = sum_choice + ((double)v1[i])/((double)n1);
			
							// If found first allele to use
							if(rand_num <= sum_choice) {
								choice1 = i;
								v1[i]--;
								n1--;
								break;
							} 
						}

						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
						sum_choice = 0.0;
						// Look or non-empty cells
						for(i = 1; i <= n; i++) {
							sum_choice = sum_choice + ((double)v1[i])/((double)n1);
			
							// If found first allele to use
							if(rand_num <= sum_choice) {
								choice2 = i;
								v1[i]--;
								break;
							} 
						}

						v1[choice1+choice2]++;
					}
					else if(rand_num <= v_c1 + v_m1) {
						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
						sum_choice = 0.0;
						// Look or non-empty cells
						for(i = 1; i <= n; i++) {
							sum_choice = sum_choice + ((double)v1[i])/((double)n1);
			
							// If found first allele to use
							if(rand_num <= sum_choice) {
								v1[i]--;
								n1--;
								v2[i]++;
								n2++;
								break;
							} 
						}
					}
					else if(rand_num <= v_c1 + v_m1 + v_c2) {
						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
						choice1 = -1;
						choice2 = -1;
						sum_choice = 0.0;
						// Look or non-empty cells
						for(i = 1; i <= n; i++) {
							sum_choice = sum_choice + ((double)v2[i])/((double)n2);
			
							// If found first allele to use
							if(rand_num <= sum_choice) {
								choice1 = i;
								v2[i]--;
								n2--;
								break;
							} 
						}

						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
						sum_choice = 0.0;
						// Look or non-empty cells
						for(i = 1; i <= n; i++) {
							sum_choice = sum_choice + ((double)v2[i])/((double)n2);
			
							// If found first allele to use
							if(rand_num <= sum_choice) {
								choice2 = i;
								v2[i]--;
								break;
							} 
						}

						v2[choice1+choice2]++;
					}
					else{//else if(rand_num <= v_c1 + v_m1 + v_c2 + v_m2) {
						rand_num = (double) rand() / RAND_MAX; // Generate uniform[0,1]
						sum_choice = 0.0;
						// Look or non-empty cells
						for(i = 1; i <= n; i++) {
							sum_choice = sum_choice + ((double)v2[i])/((double)n2);
			
							// If found first allele to use
							if(rand_num <= sum_choice) {
								v2[i]--;
								n2--;
								v1[i]++;
								n1++;
								break;
							} 
						}
					}
								
					// If simulation is over
					if(v1[n] == 1 || v2[n] == 1) {
						break;
					}	
				}
				
				// Add up the times weighted by the probability of choose k alleles for population 1
				for(i = 1; i <= n-1; i++) {
					tot_times[i] = tot_times[i] + exp(log_binom[n][k] + ((double)k)*log(x) + ((double)n-k)*log(1.0-x)) * times[i];
				}		
			}
		}

		total = 0.0;
		for(i = 1; i <= n - 1; i++) {
			tot_times[i] = tot_times[i] / ((double)num_reps);
			total = total + tot_times[i];
		}

		fprintf(outfile, "%lf", R);
		for(i = 1; i <= n - 1; i++) {
			fprintf(outfile, " %lf", tot_times[i]);
		}
		fprintf(outfile, "\n");
		printf("\tL = %lf\n", total);
		
	}
	
	fclose(outfile);

	free(tot_times);
	free(times);
	free(v1);
	free(v2);
}

double estimateDivergence(double totNumSites, double user_theta)
{
	int i = 0;
	int j = 0;
	double height = 0.0;
	int diffs = 0;
	double divergence = 0.0;
	double outgroupBranchLen = 0.0;
	double *fract = (double *)malloc((n_max + 1) * sizeof(double));
	double **derivedSpectra;
	double *thetaPi = (double *)malloc((n_max + 1) * sizeof(double));
	double temp = 0.0;

	derivedSpectra = (double **)malloc((n_max+1) * sizeof(double *));
	for(i = 0; i <= n_max; i++) {
		derivedSpectra[i] = (double *)malloc((i+1) * sizeof(double));
	}

	for(i = 0; i <= n_max; i++) {
		fract[i] = 0.0;
		thetaPi[i] = 0.0;

		for(j = 0; j <= i; j++) {
			derivedSpectra[i][j] = 0.0;
		}
	}

	// Get fraction of each sample size
	for(i = 0; i < num_sites; i++) {
		if(data[i].x == 0) {
			diffs++;
		}

		fract[data[i].n] = fract[data[i].n] + 1.0;

		derivedSpectra[data[i].n][data[i].n - data[i].x] = derivedSpectra[data[i].n][data[i].n - data[i].x] + 1.0;
	}

	for(i = 2; i <= n_max; i = i + 2) {
		fract[i] = fract[i] / ((double)(num_sites));
	}

	// Calculate Tajima's thetas
	for(i = 2; i <= n_max; i = i + 2) {
		temp = 0.0;
		for(j = 1; j < i; j++) {
			temp = temp + ((double) j * (i - j) * derivedSpectra[i][j]);
		}
		thetaPi[i] = (2.0 / ((double) i * (i - 1))) * temp;

		height = height + fract[i] * thetaPi[i] * (1.0 - 1.0/((double)i)); // Changed to include the sample size
	}

	// Delete variables from memory
	free(fract);
	free(thetaPi);

	for(i = 0; i <= n_max; i++) {
		free(derivedSpectra[i]);
	}
	free(derivedSpectra);

	// ESTIMATE OF HEIGHT = H / THETA = H / [4Nul]
	height = log(height) - log(user_theta/2.0) - log((double)totNumSites);	// Changed to divide by theta/2 rather than theta
	
	// ESTIMATE OF OUTGROUP BRANCH CONTRIBUTION = D / (THETA/2) = D / (2Nul)
	outgroupBranchLen = log((double)diffs) - log(user_theta/2.0) - log((double)totNumSites);	
	
	// ESTIMATE OF DIVERGENCE = (H + D) / 2
	divergence = 0.5*(exp(height) + exp(outgroupBranchLen));

	return divergence;
}

void printDivergence(char *outfn, double divergence)
{
	FILE *outfile = fopen(outfn, "w");

	fprintf(outfile, "%lf\n", divergence);

	fclose(outfile);
}

double *calcPolyFreq(void)
{
	int i = 0;
	double *freq = malloc((n_max+1)*sizeof(double));
	int *num = malloc((n_max+1)*sizeof(int));
	int *denom = malloc((n_max+1)*sizeof(int));

	for(i = 0; i <= n_max; i++) {
		num[i] = 0;
		denom[i] = 0;
	}

	for(i = 0; i < num_sites; i++) {
		if((data[i].x > 0) && (data[i].x < n_max)) {
			num[data[i].n]++;
		}

		denom[data[i].n]++;
	}

	for(i = 0; i <= n_max; i++) {
		if(denom[i] > 0) {
			freq[i] = ((double) num[i]) / ((double) denom[i]);
		}
		else {
			freq[i] = 0.0;
		}
	}

	return freq;
}

double *calcSubsFreq(void)
{
	int i = 0;
	double *freq = malloc((n_max+1)*sizeof(double));
	int *num = malloc((n_max+1)*sizeof(int));
	int *denom = malloc((n_max+1)*sizeof(int));

	for(i = 0; i <= n_max; i++) {
		num[i] = 0;
		denom[i] = 0;
	}

	for(i = 0; i < num_sites; i++) {
		if(data[i].x == 0) {
			num[data[i].n]++;
		}

		denom[data[i].n]++;
	}

	for(i = 0; i <= n_max; i++) {
		if(denom[i] > 0) {
			freq[i] = ((double) num[i]) / ((double) denom[i]);
		}
		else {
			freq[i] = 0.0;
		}
	}

	return freq;
}

void printFractionPolySubs(char *outfn, double* polyFreq, double *subsFreq)
{
	int i;
	FILE *outfile = fopen(outfn, "w");

	for(i = 0; i <= n_max; i++) {
		fprintf(outfile, "%d\t%lf\t%lf\n", i, polyFreq[i], subsFreq[i]);
	}


	fclose(outfile);
}

double** obtain_empirical_spectra(void)
{
	int i, j;
	double **spectra;
	int **counts;
	int sum = 0;

	spectra = (double **)malloc((n_max+1) * sizeof(double *));
	counts = (int **)malloc((n_max+1) * sizeof(int *));


	for(i = n_min; i <= n_max; i++) {
		spectra[i] = (double *)malloc((i+1) * sizeof(double));
		counts[i] = (int *)malloc((i+1) * sizeof(int));
	}

	for(i = n_min; i <= n_max; i++) {
		for(j = 1; j <= i - 1; j++) {
			counts[i][j] = 0;
			spectra[i][j] = 0.0;
		}
	}

	for(i = 0; i < num_sites; i++) {
		if((data[i].x > 0) && (data[i].x < data[i].n)) {
			counts[data[i].n][data[i].n - data[i].x]++;	// CONVERT TO DERIVED FREQUENCY SPECTRUM
		}
	}

	for(i = n_min; i <= n_max; i++) {
		sum = 0;

		for(j = 1; j <= i - 1; j++) {
			sum = sum + counts[i][j];
		}

		for(j = 1; j <= i - 1; j++) {
			if(sum > 0) {
				spectra[i][j] = counts[i][j] / ((double) sum);			
			}
		}
	}


	for(i = n_min; i <= n_max; i++) {
		free(counts[i]);
	}

	free(counts);

	return spectra;
}

void write_empirical_spectra(char *outfn, double **spectra)
{
	FILE* outfile;
	int i, j;

	outfile=fopen(outfn, "w");

	fprintf(outfile, "%d %d\n", n_min, n_max);
	for(i = n_min; i <= n_max; i++) {
		for(j = 1; j <= i - 1; j++) {
			fprintf(outfile, "%lf ", spectra[i][j]);
		}
		fprintf(outfile, "\n");
	}


	fclose(outfile);
}

void free_empirical_spectra(double **spectra)
{
	int i;

	for(i = n_min; i <= n_max; i++) {
		free(spectra[i]);
	}

	free(spectra);
}


void usage() 
{
	printf("\n\n*****************************************\n");
	printf("***************** USAGE *****************\n");
	printf("*****************************************\n\n");
	
	printf("********** Prepare scan **********\n");
	printf("Get mean interspecies coalescence time: ./BALLET -inter_coal_time CombinedSNPFile TotNumSites 4Nu DivFile\n");
	printf("Get proportion of polymorphic and subsituted sites: ./BALLET -poly_sub CombinedSNPFile PolySubFile\n");
	printf("Get empirical spectra: ./BALLET -spect CombinedSNPFile SpectFile\n");
	
	printf("\n********** Balancing selection **********\n");
	//printf("Generate spectrum: ./BALLET -SimSpect n x theta1 theta2 num_replicates OutFile\n");
	printf("Scan with T1: ./BALLET -T1 WINDOWSIZE DivFile PolySubFile SNPFile RecFile OutFile\n");
	printf("Scan with T2: ./BALLET -T2 WINDOWSIZE DivFile PolySubFile SpectFile SNPFile RecFile PATH OutFile\n");
	
	exit(-1);
}
