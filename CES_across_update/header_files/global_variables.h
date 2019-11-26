
// Declaring Header Files
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>//for mkdir
#include <sys/types.h>//for types w/in mkdir
#include <math.h>
#include <time.h>
#include <sys/time.h> //gettimeofday()
#include <string.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>//uniform rng
#include <gsl/gsl_randist.h>//rng from distns


//Preprocessor directives
#define MAX_DDIM 30000
#define MAX_AA 6500 
#define MAX_NTS 19500 
#define MAX_LOCI 6000 //set max loci to read in and simulate
#define MAX_TIME -20//Max number of evolutionary steps to simulate, 
// if set to a negative value it will calculate the appropriate max time to have X substitutions/nt.
//no longer used.  Defined via coommand line
#define NE 1.36e7 //effective population size.
//#define Q 1.519e-5 //fitness scaling coefficient
#define MU 1E-9 //per generation mutation rate

//Turn off error codes
int gsl_errors_off = 0;
int print_debug = 0;
 
// Declaring Global Variables
double 	A1; //Cost of initiation of translation in units of ATP. Default val set in command line to  A1 = 2;
double	A2; //Cost each elongation step in ATPs, default set in command line to A2=4
double  B; // Nonsense error rate. units of 1/sec. 
double	pi=M_PI; //
double mu = MU; //per nt mutation rate.
double  min_z_max = 1e-10; //minimum value to start integration routines with;
double max_z_max = 1e100; //set an upper bound for z
double z_max_factor = 10; //amount to multiply z_max by when searching for value
double relerr = 1e-6; //relative error for integration routines
double abserr = 1e-20; //absolute error for integration routines.  aka tolerance
double at_bias = 0.5; //parameter for adjusting nt composition due to biased mutation and/or gene conversion.
                      //According to table 6.1 in Lynch (2007) on p. 125, the observed AT bias is 0.62
double mu_bias = 1; //parameter calculated based on at_bias.  mu_bias = at_bias/(1-at_bias)
double gamma_ratio = 1; //ratio of transition mutation rate to transversion mutation rate 
                        //This is \alpha/\beta in the Tamura and Nei (1993) where we assume \alpha_1 = \alpha_2
                        //Calculations using data from Lynch give value of approximately 1.22 for S.c. 
double mutation_matrix[4][4]; //mutation matrix order is ATCG.  See mutation.matrices.nb for details.

double global_max_time=MAX_TIME; //Max number of evolutionary steps to simulate, if set to a negative value it will set the max time so that, under neutrality, there are on average X substitutions/nt.
double Ne = NE;
double Q = 4.19e-7; // Fitness scaling coefficient

int ignore_aa;  //decrease true aa_count by this amount.
int random_start=1; //indicate whether or not to start with a random sequence (1) or the seq read in (0);
int n_aa=0,//Number of aa in a sequence
  n_seq, //set when reading fasta file
  pout, //print out configuration flag.  Stores -W arguement passed on command line
  pconf,
  gconf,
//  calcamean=0,//indicates whether to calculate arithmetic mean.  Will be reset based on arguments passed
//  calcgmean=0,
//  calcvar=0,
//  calcmode = 1,//indicates which method to use to estimate the posterior mode. 1:gamma based 2: numerical search
  runs,
  analytic,
  max_aa,
  MLE;
int print_evol_fasta = 0; //indicate whether or not to print seq evolution to fasta file
int print_delta_eta_vals = 0; //indicate whether or not to print list of delta eta vals for each step
//can be set to 1 via command line

//default tRNA file
char def_tRNA[150]="tRNA_files/S.cerevisiae.tRNA.tsv";
char fasta_file[150] ="./fasta/S.cerevisiae.S288c.with.phi.fasta";
char out_prefix[150], trna[150];
char out_folder[150];//derived from out_prefix
char command_line[1000];//save command line for printing
char exec_time[360];//save time command executed at.


//  time_t start,comp_stop, print_stop;
time_t start;
struct tm * timeinfo;
struct timeval time_val, comp_start, comp_stop, print_stop, time_diff;

// Defining structure for a single amino acid
// keep track of E(1/c), E(1/c^2), E((c+b)/c)... 
// which are expectations wrt the translation rates for the aa's set of codons
struct amino_acid
{	
  char codon[6][4];
  char aa;
  int num_codons; //was cc
  int codon_index[6]; //used when generating random sequences
  //expected values for elongation related rates
  double e_invc; // 1/c
  double e_cpbinvc;//(c+b)/c
  double e_invc2;//1/c^2
  double e_cpb2invc2;//(c+b)^2/c^2
  double e_cpbinvc2;//(c+b)/c^2 --used in calc cov b/w codons
  double elong_rate[6]; //elongation rate of codons for AA, was tr_rate
  double elong_pr[6];//Pr of successful elongation c/(c+b)
  double fail_pr[6];//Pr of unsuccessful elongation b/(c+b)
  double max_elong_rate;//max elongation rate for the amino acid
  double max_elong_pr;//max elongation rate for the amino acid
  double min_elong_rate;//min elongation rate for the amino acid
  double min_elong_pr;//min Pr of successful elongation for the amino acid
  double neutral_obs_pr[6];//Pr of observing this codon based solely on its nt sequence and any AT bias.
  double neutral_obs_pr_cum[6];//Cumulative probability of observing this codon. Useful for when choosing a random codon when there is AT bias
};



struct codon_struct
{	
  char codon[4];
  char aa;//aa letter
  
  unsigned short int synonym_index[5];//Alternative synonyms for the same aa listed by their codon index #
  unsigned short int one_step_synonym_index[5];//as above but subset that is one step mutant from codon.
  
  unsigned short int aa_index;
  
  unsigned short int num_synonym;
  unsigned short int num_one_step_synonym;

  
  double  elong_pr_ratio[6];//(c_i/(c_i+b))/(c_j/(c_j+b)) for one step mutant from codon c_i.
  double  delta_b_over_c[6];//(b/c_i - b/c_j) for one step mutant from codon c_i.
  double  one_step_relative_mutation_rate[5];//mutation rate to one step neighbors of current codon
  //  double elong_rate; //elongation rate of codons for AA, was tr_rate
  //as above but subset that is one step mutant from codon.

  double b_over_c;//b/c_i
  double elong_rate; //elongation rate of codons for AA, was tr_rate
  double elong_pr;//Pr of successful elongation c/(c+b)
  double fail_pr;//Pr of unsuccessful elongation b/(c+b)
};


struct codon_struct Codon[64];


struct seq_struct //it might be better to separate sequence info from stats from SEMPPR
//by  creating a separate structure class
{	
  char id[16]; //Gene ID or ORF name
  char codon_seq[MAX_AA][4];//codon sequence
  //save memory by using unsigned short ints
  unsigned short int codon_index[MAX_AA]; //codon index code
  unsigned short int aa_index[MAX_AA]; //aa index code
  
  int aa_count; //aa count
  double sigma_vec[MAX_AA];//vector of \sigma(i) values.  
  //Note that sigma_0 repesents the start codon so its value is 1. 
  double sigma_ratio_vec[MAX_AA];//vector of \sigma(i)/\sigma(n) values
  //  double xi_vec[MAX_AA];//vector of 1/(1-\sigma(n)) \sum_{i=j}^n (a1+a2 i) \sigma_{n-1} b/(c_i+b) values
  //i.e. xi_vec[0] = \xi, xi_vec[1] = \xi - (a1 +a2)(b/(c_1+b)) or \sum_{i=2}^n (a1+a2 i) \sigma_{n-1} b/(c_i+b) values
  double max_time;
  int ARR;
  double pA; //alpha parameter for f(eta)~beta distribution
  //used by gsl
  double pB; //beta+alpha parameter for f(eta)~beta distribution
  double pD; //eta_diff = eta_max - eta_min
  double pN; //eta_obs-eta_min
  double pdenom; //Scaling term for pdf
  double MEAN; //arithmetic mean expression ?
  double z_max;
  double z_mode;//mode of z based on eta~beta dist
  double z_prcntl[202]; 
  double z_amean;
  double z_gmean;
  double z_var;
  double z_sd;
  double z_mode_gamma; //mode based on assuming eta ~ gamma
  double sigma_obs; //obs sigma value
  double xi_obs; //obs xi value
  double eta_obs; //obs eta value
  double eta_min; //min eta value
  double eta_max;//max eta value
  double eta_mean; //mean of eta distn
  double eta_var; //var of eta distn
  double eta_sd; //sd of eta distn

  double alpha;//alpha parameter for f(eta)~beta distribution
  double beta;//beta parameter for f(eta)~beta distribution
  double gamma_alpha;//alpha parameter for f(eta)~gamma distribution
  double gamma_beta;//beta parameter for f(eta)~gamma distribution

  // following added for evolution simulations
  double eta_initial;//eta value at start of simulation
  double phi_obs; //assumed production rate
  double z_obs; //assumed production rate, scaled by Ne and q
  double delta_eta_mean;//mean of distn of delta_eta values for wt allele
  double delta_eta_var; //var of distn of delta_eta values for wt allele

  double evol_delta_eta_mean;//mean diff b/w wt_eta and mut_eta;
  double evol_delta_eta_var; //var b/w wt_eta and mut_eta;
  double evol_time; //var b/w wt_eta and mut_eta;
  
  int evol_steps; //number of times wt allele is replaced by a mutant
  int num_synon_mut; //num of mut that were synonymous


  //  double z_lowcut;
  //double z_upcut;
};

//uint8_t  is an unsigned, 8 bit type that ranges from 0-255;
//an alternative form is unsigned char 
//unsigned short int should range from 0 to 65535

struct seq_history_struct //structure to record the evolutionary history of a sequence
{	
  char id[16]; //Gene ID or ORF name
  char wt_codon_seq[MAX_AA][4];//codon sequence
  //  unsigned short int  mut_codon_pos[MAX_EVOL_STEPS]; //record codon positions where mutations occur
  //unsigned short int  mut_nt_pos[MAX_EVOL_STEPS]; //record nt positions w/in codon where mutations occur
  //unsigned short int  mut_codon_index[MAX_EVOL_STEPS];//list of codons mutated. 
  //unsigned short int synon_mut[MAX_EVOL_STEPS];//indicates whether mutation is synonymous (1) or not (0)
  //unsigned short int fixation[MAX_EVOL_STEPS];//indicates whether fixation occurs (1) or not (0)
  
  //double sigma_n_vals[MAX_EVOL_STEPS];
  //double xi_vals[MAX_EVOL_STEPS];
  //double eta_vals[MAX_EVOL_STEPS];

  double eta_min;
  double eta_max;
  double eta_mean;
  double eta_var;
};


struct 	amino_acid  AA[22];
struct 	seq_struct Seq[MAX_LOCI];
