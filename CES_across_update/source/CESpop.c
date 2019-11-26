/* Software Description


   Codon Evolution Simulation (CES)
   Stochastic simulation of the evolution of a codon sequence under selection against nonsense errors.

   This software was developed as part of

	"Gilchrist, M.A., P. Shah, and R. Zaretzki (2009 or 2010 depending on final publication date). 
	Measuring and detecting molecular adaptation in codon usage against nonsense errors 
	during protein translation. Genetics"

This work should be cited when this software is modified or used.


*/

/* Software License

   CES is free software: you can redistribute it and/or modify it under the terms of the 
   GNU General Public License (GPL) Version 2 as published by the Free Software Foundation 
   at http://www.gnu.org/licenses.

   This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
   for more details. 

   Users are free to distribute and modify this software according to this license so long as this and all of the
   statements above  are not removed or modified.

*/

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
#define PHI 1.0 //scales phi by this value

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
double Phi = PHI;
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

// Error message when improper use of the code
int wrong()
{	printf(	"\nError in command line:\nCES Exiting...\n");
  exit(1);
}




// Opening a tRNA file
int Read_tRNA_File(char *filename)//, int *num_aa, int *num_codons)
{	
  int i, j, k, jmax, imax;
  FILE *file_handle;
  char curr_char;
  char aa_list[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'Z', 'S', 'T', 'V', 'W', 'Y', 'X'}; //'X' represents the 'stop' amino acid, 'Z' represents the smaller set of serine codons: AGT and AGC.  The other ser codons are TCT, TCC, TCG, and TCA.
  //Note program expects 'X' to be at end of list so 'Z' is put in internally.
  double elong_rate;
  int num_codons;
  int aa_processed;
  int tmp;
  int max_aa = 22; //20 AA but serine may get split and we may have stop codons


  j=0;
  //initialize AA structures
  for(i=0; i<max_aa; i++){
    AA[i].aa = aa_list[i];
    AA[i].num_codons = 0;
  }

  file_handle=fopen(filename, "r");
  if(!file_handle)
    {
      printf("\ntRNA File: %s Doesn't Exist\ntRNA", filename);
      wrong();
    }
		
  //get aa index 
  curr_char=fgetc(file_handle);
  
  while(curr_char != EOF)
    {	
      i = 0; //set AA index to 0
      //flip through AA until you get a match
      while(curr_char != AA[i].aa && i < max_aa)
	{
	  i++;
	}

      if(i ==max_aa) {
	printf("\nAA index in file did not match any known AA. Exiting...\n");
	exit(1);
      }
      else{
	//get current codon count
	num_codons = AA[i].num_codons;

	if(num_codons ==6) //check num_codons value 
	  {
	    printf("\nCodon count for AA is greater than maximum possible value of 6.  Exiting...\n");
	    exit(1);
	  }
      }	
      
      curr_char=fgetc(file_handle);//get next character. Should be a \t

      //read in codon
      for(j=0;j<3;j++)//load codon sequence and put in codon_index codon index
	{	
	  AA[i].codon[num_codons][j]=fgetc(file_handle);
	}
      AA[i].codon[num_codons][j]='\0';


      curr_char=fgetc(file_handle);//get next char. Should be a \t

      //read in elongation rate
      tmp = fscanf(file_handle,"%lf",&AA[i].elong_rate[num_codons]);
      
 
      //elong_rate = AA[i].elong_rate[num_codons];
      //increment codon count for the aa
      AA[i].num_codons++;
      num_codons++;


      //read until end of line or EOF
      while((curr_char !='\n')&&(curr_char !=EOF)) curr_char=fgetc(file_handle);

      //get aa index for next codon if not at EOF
      if(curr_char !=EOF) curr_char=fgetc(file_handle);

    }
  fclose(file_handle);


  //check to make sure the correct # of AA and codons are defined
  aa_processed=0;
  for(i=0; i<max_aa;i++)
    {
      if(AA[i].num_codons>0){
	aa_processed++;
      }
    }
  return aa_processed;//return # of 

}



int Read_FASTA_File(char *filename)
{	

  int i,j,k,h,l,n,r,id;
  char curr_char, codon[4];
  char curr_line[82];
  char *char_ptr;
  char *str_test;
  char nt_seq[MAX_NTS];
  int end_of_seq;
  int imax, jmax, kmax, lmax,nt_count, aa_count;
  FILE *fh;
	
  fh=fopen(filename, "r");
  if(!fh)
    {	
      printf("\nFASTA/Sequence File Doesn't Exist\n");
      wrong();
    }
  codon[3]='\0';

  //  curr_char='';
  str_test=curr_line;

  //Read until you get to a new line with a >
  do{
    str_test = fgets(curr_line, 82,fh);
  }while(curr_line[0]!='>' && str_test!=NULL);

  i=0;
    
  while(str_test!=NULL && i < MAX_LOCI)
    {

      //line should be start of new sequence
      //read in seq id
      k=0;
      while((curr_char!=' ')&&(curr_char!='\t')&&(curr_char!='\n')&&(k<15))
	{
	  curr_char = curr_line[k+1];//offset by 1 b/c of >
	  Seq[i].id[k]=curr_char;
	  k++;
	}
      if(k==15){
	printf("Seq ID greater than 15 characters. Exiting");
	exit(1);
      }      
      //replace new line with null character
      Seq[i].id[--k]='\0';

      //Read in phi value
      //locate start of text  "phi = "
      char_ptr = strstr(curr_line, "phi=");//removed space b/w phi and = sign 

      if(char_ptr==NULL){
	printf("Phi value not found for seq %d (%s). Exiting.\n", i, Seq[i].id);
	exit(1);
      }
      
      char_ptr+=strlen("phi =");//go to end of match string
      //note strtod will skip any whitespace preceeding the number
      Seq[i].phi_obs = Phi*strtod(char_ptr, NULL);//first argument is ptr to str location of double
      // arguementsecond is a pointer to the place afterwards

      
      nt_count = 0;
      nt_seq[0]='\0';

      
      end_of_seq=0;

      
      //get next line
      str_test = fgets(curr_line, 82,fh);

      while(str_test!=NULL && end_of_seq == 0){
	//test to make sure string isn't too long
	if(strlen(curr_line)>81){
	  printf("ERROR: FASTA file column is >80 characters. Exiting");
	  exit(1);
	}
     
	//read through until you hit a '>'
	curr_char = curr_line[0];

	switch(curr_char){
	case 'A':
	case 'T':
	case 'G':
	case 'C':
	case 'a':
	case 't':
	case 'g':
	case 'c':
	  k=strlen(curr_line);
	  if(curr_line[k-1]=='\n'||curr_line[k-1]==EOF) k--;
	  nt_count+=k;
	  if(nt_count >= MAX_NTS){
	    printf("Seq[%d].id = %s has %d nts.  MAX_NTS is %d. Exiting", i, Seq[i].id, nt_count, MAX_NTS);
	    exit(1);
	  }
	  strncat(nt_seq, curr_line, k);
	  //get next line
	  str_test = fgets(curr_line, 82,fh);
	  break;
	case '>':
	  end_of_seq=1;
	  //don't get next line
	  break;
	case '\n':
	  //get next line
	  str_test = fgets(curr_line, 82,fh);
	  break;
	case EOF: //str_test!=NULL should prevent reaching this point
	  end_of_seq=1;
	  break;
	default:
	  printf("Read in unexpected character (%c) at beginning of line. Exiting.", curr_char); 
	  exit(1);
	  break;
	}
      }

      //process nt seq
      if(nt_count % 3 != 0){
	printf("Seq[%d].id = %s had %d nts, which is not a multiple of 3. Exiting", i, Seq[i].id, nt_count);
	exit(1);
      }
      else{//everything seems okay
	aa_count = (nt_count/3)-1;//subtract 1 b/c of stop codon
	Seq[i].aa_count=aa_count;
	
	k=0;
	
	for(j=0;j<=aa_count;j++){//include stop codon
	  Seq[i].codon_seq[j][0]=nt_seq[k++];
	  Seq[i].codon_seq[j][1]=nt_seq[k++];
	  Seq[i].codon_seq[j][2]=nt_seq[k++];
	  Seq[i].codon_seq[j][3]='\0';
	}
	
      }
    
      i++;//increment seq id
	
    }
  return i;
}





void Convert_Codon_Seq_to_Codon_Index(char codon_seq[][4], unsigned short int *index_vec, int aa_count)
{
  int i, j, k, match;
  //  char codon[4];

  for(i=0; i<aa_count; i++){
    
    //    strcpy(codon, codon_seq[i]);

    j=0;
    do{
	match = strcmp(Codon[j++].codon, codon_seq[i]);
    }while(match!=0 && j<64);

    if(j==64){
      printf("ERROR: Can't match %s to a codon. Exiting", codon_seq[i]);
      exit(1);
    }
    else{
      index_vec[i]=--j;//counter act j++ in do loop
    }
    
  }
}

int Calc_Dimensionality(unsigned short int *codon_index_vec, int aa_count){
  int i, j, k, Dim;
  
  Dim=0;
  for(i = 1;i<aa_count;i++){
    Dim+=Codon[codon_index_vec[i]].num_one_step_synonym;
  }
  
  return Dim;


}

int Generate_D_Array(unsigned short int *codon_index_vec, unsigned short int *D_array[][2], int aa_count){
  int i, j, k, jmax, imax;
  int D_index;
  short unsigned int codon_index;


  //  imax = (unsigned short int)(aa_count);
  

  D_index=0;
  for(i = 1;i<aa_count;i++){
    codon_index =codon_index_vec[i];
    jmax = Codon[codon_index].num_one_step_synonym;

    //go through each synonym
    for(j=0;j<jmax;j++){
      *(D_array[D_index][0])= (unsigned short int)(i);
      *(D_array[D_index++][1])= (Codon[codon_index].one_step_synonym_index[j]);
    }
    
  }
}

void Generate_Random_Codon_Seq_and_Index(char codon_seq[][4], unsigned short int *codon_index, unsigned short int *aa_index, int aa_count){
  int i, j, k, h, num_codons, test;
  double rnd;

  //things related to rng
  struct tm * timeinfo;
  struct timeval curr_time;

  const gsl_rng_type * rng_type;
  gsl_rng * rn;  


  //Set up GSL RNG
  gsl_rng_env_setup();
  

  
  //Create a generator chosen by the 
  //environment variable GSL_RNG_TYPE
  //use default gsl for generating uniform rn 
  //from which other rn functions are derived 

  rng_type = gsl_rng_default;
  rn = gsl_rng_alloc (rng_type);

  //seed rng using clock
  gettimeofday(&curr_time,NULL);
  if(print_debug){
    printf("Seeding RNG with clock. Current time is %d %d\n", (int)(curr_time.tv_sec), (int)(curr_time.tv_usec));
    fflush(stdout);
  }

  gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));// (const gsl_rng * r, unsigned long int s)   

  //end RNG set up

  
  for(i=1;i<aa_count;i++)
    {  
      k=aa_index[i];
      num_codons = AA[k].num_codons;
      if(num_codons > 1){
	rnd= gsl_rng_uniform(rn);
	test =0;
	for(j=0;j<num_codons && test==0;j++){  
	  if(rnd<AA[k].neutral_obs_pr_cum[j]) //Picks a codon based on its neutral probability and takes into account AT bias
	    {  
	      //copy codon index and seq
	      codon_index[i]=AA[k].codon_index[j];
	      strcpy(codon_seq[i],AA[k].codon[j]);
	      test=1; //break;
	    }
	}
      }
    }
}


void Generate_Mutation_Matrix(double mut_matrix[][4], double at, double g_ratio)
{
  int i, j, imax = 4;

  //define based on AT bias
  for(i=0; i<2; i++){
    for(j=0; j< imax; j++){
      mut_matrix[i][j] = at;
    }
  }  
  for(i=2; i<imax; i++){
    for(j=0; j< imax; j++){
      mut_matrix[i][j] = 1-at;
    }
  }

  //Adjust for transition vs. transversion
  //note these are on the off diagonal
  j = 3;
  for(i=0; i< imax; i++){
    mut_matrix[i][j] *= g_ratio;
    j--;
  }

  
  //For completeness set diagonal values to zero
  for(i=0; i< imax; i++){
    mut_matrix[i][i] =0;
  }

}


int Calc_Factor_Vec(unsigned short int *codon_index_vec, double *factor_vec, int aa_count){
  int i, j, k, D;
  double factor;
  for(i = 0;i<aa_count;i++){
    Codon[codon_index_vec[i]].num_one_step_synonym;
  }
  
}


void Convert_Codon_Index_to_Sigma_Vec(unsigned short int *index_vec, double *sigma_vec, int aa_count)
{
  int i, j, k, match;
  double sigma_i;
  char codon[4];

  sigma_vec[0]=1.0;
  sigma_i=1.0;

  for(i=1; i<aa_count; i++){//skip first amino acid
    
    sigma_i*=Codon[index_vec[i]].elong_pr;
    sigma_vec[i] = sigma_i;
  }

}


void Calc_B_over_C_Vec(unsigned short int *codon_index_vec, double *b_over_c_vec, int aa_count){
  int i, j,k;

  for(i=0;i<aa_count;i++){
    b_over_c_vec[i] = B/Codon[codon_index_vec[i]].elong_rate;
  }

}




void Convert_Sigma_Vec_to_Sigma_Ratio_Vec(double *sigma_vec, double *sigma_ratio_vec, int aa_count)
{
  int i, j, k, match;
  double sigma_n;
  char codon[4];

  if(print_debug){
    printf("Entering Convert_Sigma_Vec_to_Sigma_Ratio_Vec");
    fflush(stdout);
  }
    sigma_n=sigma_vec[aa_count-1];
  for(i=0; i< aa_count; i++){//don't think the first is every used
   sigma_ratio_vec[i] = sigma_vec[i]/sigma_n;
  }

}

//recalcs sigma vec given that elong_pr of codon at position i changes by a factor
void Recalc_Sigma_Vec(double *sigma_vec, int mut_codon_pos, double factor, int aa_count)
{
  int i, j, k, match;
  char codon[4];


  for(i=mut_codon_pos; i< aa_count; i++){
    sigma_vec[i] *= factor;
  }

}

void Convert_Codon_Index_to_AA_Index(unsigned short int *codon_index_vec, unsigned short int *aa_index_vec, int aa_count)
{
  int i, j, k, match;
  for(i=0; i< aa_count; i++){
    aa_index_vec[i]=Codon[codon_index_vec[i]].aa_index;
  }

}

// Calculating Xi based on a vector of sigma values
double Calc_Xi(double *sigma_vec, int aa_count)
{      
  int i=0;
  double xi=0.;
  double sigma_n;
  
  sigma_n = *(sigma_vec+(aa_count-1));
  sigma_vec++;

  //xi=(A1+A2)*(1-*sigma_vec);//should be zero
  for(i=1;i<aa_count;i++)
    {	
      
      xi+=(A1+A2*(i))*(*(sigma_vec-1)-*(sigma_vec++));
      //\sigma_{i-1} b/(b+c_i) = \sigma_{i-1}-\sigma_i 
      //modified from A2*(i+1)
      //rationale: should be i-1, but C indexing offsets things by 1, so should just be i
    }

  xi *=(1/(1-sigma_n));

  return xi;
}


//Calc_Eta 
double Calc_Eta(double sigma_n, double xi, int aa_count) 
{
  double eta;
  eta=(1- sigma_n)*xi/sigma_n+(A1+A2*aa_count);
  return eta;
}


void Calc_Delta_Eta_First_and_Second_Term_Vecs(double *sigma_ratio_vec, double *b_over_c_vec, double *delta_eta_first_term_vec, double *delta_eta_second_term_vec, int aa_count){
  int i, j, k;
  double curr_cost;
  double sigma_n, delta_eta_first_term, delta_eta_second_term;
  //  double z; 


  curr_cost = (double)(A1);

  *(delta_eta_first_term_vec++)=0.0;//zero b/c there is no chance of a nonsense error at the first codon.
  *(delta_eta_second_term_vec++)=0.0;//zero b/c there is no chance of a nonsense error at the first codon.

  b_over_c_vec++;
  sigma_ratio_vec++;


  delta_eta_first_term = 0.;

  for(i=1;i<aa_count;i++){
    //term is a summation from j=0 to i-1, so we don't update the first term until the end
    *(delta_eta_first_term_vec++)=delta_eta_first_term;

    curr_cost+=A2;

    *(delta_eta_second_term_vec) = curr_cost*(*(sigma_ratio_vec++));

    delta_eta_first_term+=(*(b_over_c_vec++))*(*(delta_eta_second_term_vec++));

    //(A1 + A2 i)*b_over_c_vec[i]*sigma_ratio_vec[i];
  }

  
}

void Calc_Delta_Eta_Vec(unsigned short int *codon_index_vec, double *b_over_c_vec, double *first_term_vec, double *second_term_vec, double *delta_eta_vec, int aa_count, int D_dim){
  int i, j, k, jmax, D_index;
  double curr_cost;
  unsigned short int codon_index, synon_index;
  double delta_eta, first_term, second_term;
  double z;  //z = \frac{\?igma_{n,i}}{\sigma_{n,j}} =  \left(\frac{c_{k,i}}{c_{k,i}+b} \frac{c_{k,j}+b}{c_{k,j}}\right)

  D_index=0;

  delta_eta = 0;

  curr_cost = (double)(A1);
  codon_index_vec++;
  first_term_vec++;
  second_term_vec++;
  
  //go through each codon
  for(i = 1;i<aa_count;i++){
    codon_index =*(codon_index_vec++);
    jmax = Codon[codon_index].num_one_step_synonym;
    
    curr_cost+=A2;

    first_term = *(first_term_vec++);
    second_term = *(second_term_vec++);
    //go through each synonym
    for(j=0;j<jmax;j++){
      z = Codon[codon_index].elong_pr_ratio[j];
      
      delta_eta_vec[D_index++]= first_term *(1-z) + second_term *Codon[codon_index].delta_b_over_c[j]; 
    }
    
  }
  if(D_index != D_dim){//not sure why I need to add one before comparing to D_dim
    printf("ERROR in Calc_Delta_Eta_Vec: %d out of %d indices cycled through. Exiting", D_index, D_dim);
    exit(1);
  };
}



//Calculates pi(i->j) value according to Sella and Hirsh
void Calc_Pi_Vec(double *delta_eta_vec, double *pi_vec, double *ptr_pi_total, double phi, double qPhi, double two_qNePhi, int aa_count, int D_dim){
  int i, j, k;
  double pi, delta_eta;
  double numerator, denominator;
  double pi_total;
  double invNe;


  //assuming diploid.  If haploid need to multiply qPhi by 2
  invNe = 1./Ne;

  pi_total = 0;
  for(i=0;i<D_dim;i++){
    delta_eta = delta_eta_vec[i];
    
    if(delta_eta ==0 || phi==0){//was producing strange behavior when phi=0.  First codon was always changing.
      //pure drift process
      pi = pi_vec[i]=invNe;
       }
    else{
    //numerator = -expm1(minuxtwoqphi*delta_eta);//expm1(x) = exp(x)-1, -expm1 = (1-exp(x))
    //denominator = -expm1(minustwoqNephi*delta_eta);
    pi = pi_vec[i] =expm1(- qPhi*delta_eta)/expm1(-two_qNePhi*delta_eta); //numerator/denominator;
    } 
    pi_total+=pi;
  }
  *(ptr_pi_total) = pi_total;
}


//Weights pi by mu to get replacment pr
void Calc_Replacement_Pr_Vec(double *mutation_vec, double *pi_vec, double *pr_vec, double *ptr_pr_total, int D_dim){
  int i, j, k;
  double pr;
  double pr_total;
  double invNe;

  pr_total = 0;
  for(i=0;i<D_dim;i++){
    pr = pr_vec[i] =mutation_vec[i]*pi_vec[i];
    pr_total+=pr;
    } 

  *(ptr_pr_total) = pr_total;
}

//Find max/min elongation rates, calc elong_pr, fail_pr.
int Process_AA_Information(){
  char nt;
  char codon[4];
  int i, j, k;
  int num_codons;
  int max_aa=22;
  double sum_neutral_obs_pr;
  double elong_rate, min_elong_rate, max_elong_rate;

  mu_bias = at_bias/(1-at_bias);  //calculate biased mutation rate based on AT bias
                                  //Added 3/2/09
  //check to make sure the correct # of AA and codons are defined
  j=0;
  k=0;

  if(print_debug){
    printf("Processing AA information\n");
    fflush(stdout);
  }

  //check for if stop codons defined. If not, define them.
  if(AA[max_aa-1].num_codons ==0)
    {
    AA[max_aa-1].num_codons = 3;
    strcpy(AA[max_aa-1].codon[0],"TGA");
    strcpy(AA[max_aa-1].codon[1],"TAG");
    strcpy(AA[max_aa-1].codon[2],"TAA");
    }

  for(i=0; i<max_aa;i++)
    {
	k+=AA[i].num_codons;
      }
  
  if(k!=64){
    printf("Only %d codons defined. Expecting 64. Exiting", k);
    exit(1);
  }
 

  for(i=0; i<max_aa-1;i++)
    {
      num_codons = AA[i].num_codons;
      
      sum_neutral_obs_pr = 0; //added line 11/21/08

      for(j=0;j<num_codons;j++){
	elong_rate = AA[i].elong_rate[j];
	AA[i].elong_pr[j]=elong_rate/(elong_rate+B);
	if(print_debug){
	  printf("%d:%d %s: %f\n", i, j, AA[i].codon[j], AA[i].elong_pr[j]);
	  fflush(stdout);
	}
	AA[i].fail_pr[j]=B/(elong_rate+B);
      //first part of calculating the pr of obs codon under neutral model but with AT bias
	//set to 1 and then scale based on at_bias
	AA[i].neutral_obs_pr[j] = 1;

	strcpy(codon, AA[i].codon[j]);
	for(k=0;k<3;k++){
	  nt = codon[k];
	  if(nt=='A' || nt=='T'){
	    AA[i].neutral_obs_pr[j] *= mu_bias; // old, wrong code: at_bias;
	      }//else{
	  //AA[i].neutral_obs_pr[j] *= (1-at_bias);
	  //}
	}
	//keep track of total for an AA.  This will be used to scale previous calculations
	sum_neutral_obs_pr += AA[i].neutral_obs_pr[j];
	AA[i].neutral_obs_pr_cum[j]=sum_neutral_obs_pr;
      }

      //scale each codon by sum_neutral_obs_pr to ensure they sum to 1
      for(j=0; j<num_codons; j++){
	AA[i].neutral_obs_pr[j]/=sum_neutral_obs_pr;
	AA[i].neutral_obs_pr_cum[j]/=sum_neutral_obs_pr;
      }
    }


  if(print_debug){
    printf("\tFinding min/max elongation rates and pr\n");
    fflush(stdout);
  }
  

  
  //Assign min/max_elong_rate/pr
  for(i=0; i<max_aa-1;i++)
    {
      num_codons = AA[i].num_codons;
      
      min_elong_rate=AA[i].elong_rate[0];
      max_elong_rate=min_elong_rate;
      
      
      //start at 1 b/c already processed first codon
      for(j=1;j<num_codons;j++){
	    
	if(min_elong_rate > AA[i].elong_rate[j]) {
	  min_elong_rate = AA[i].elong_rate[j];
	}

	if(max_elong_rate < AA[i].elong_rate[j]) {
	  max_elong_rate = AA[i].elong_rate[j];
	}
	    
      }

      AA[i].min_elong_rate=min_elong_rate;
      AA[i].max_elong_rate=max_elong_rate;

      AA[i].min_elong_pr=min_elong_rate/(min_elong_rate+B);
      AA[i].max_elong_pr=max_elong_rate/(max_elong_rate+B);


    }


 

// Calculating the first and second moments for each amino acid
  if(print_debug){
    printf("calculating moments\n");
    fflush(stdout);
  }
  for(i=0;i<max_aa-1;i++){
      if(print_debug){
	printf("calculating values for aa %d\n", i);
	fflush(stdout);
      }
      AA[i].e_invc=0;
      AA[i].e_cpbinvc=0;
      AA[i].e_invc2=0;
      AA[i].e_cpb2invc2=0;
      AA[i].e_cpbinvc2=0;
      for(j=0;j<AA[i].num_codons;j++){
          AA[i].e_invc+=1/(AA[i].elong_rate[j])*AA[i].neutral_obs_pr[j];
	  AA[i].e_cpbinvc+=(B+AA[i].elong_rate[j])/(AA[i].elong_rate[j])*AA[i].neutral_obs_pr[j];
	  AA[i].e_invc2+=1/(pow(AA[i].elong_rate[j],2))*AA[i].neutral_obs_pr[j];
	  AA[i].e_cpb2invc2+=(pow(B+AA[i].elong_rate[j],2))/(pow(AA[i].elong_rate[j],2))*AA[i].neutral_obs_pr[j];
	  AA[i].e_cpbinvc2+=(B+AA[i].elong_rate[j])/(pow(AA[i].elong_rate[j],2))*AA[i].neutral_obs_pr[j];
      }
  }
    
  if(print_debug){
    printf("exiting moments\n");
    fflush(stdout);
  }

return(1); //exited successfully
}


void Generate_Codon_Structures()
{
  int num_syn, num_one_step_syn;
  int i, j, k, test;
  int num_matches;
  short unsigned int codon_index, synon_index, max_aa=22;
  int codon_match[3];//use to determine which position and nts involved in differences b/w codons.
  double sum_neutral_obs_pr;
  double relative_mut_rate;

  char codon[4];
  char synonym[4];
  char orig_nt, mut_nt;

  codon_index = 0;

  for(i=0;i<max_aa;i++){
    for(j=0;j<AA[i].num_codons;j++){
      AA[i].codon_index[j]=codon_index;
      strcpy(Codon[codon_index].codon, AA[i].codon[j]);
      Codon[codon_index].aa = AA[i].aa;
      Codon[codon_index].aa_index = i;
      Codon[codon_index].elong_rate = AA[i].elong_rate[j];
      Codon[codon_index].elong_pr = AA[i].elong_pr[j];
      Codon[codon_index++].fail_pr = AA[i].fail_pr[j];
	}
  }

  //Define synonyms
  for(i=0;i<64;i++){
    num_syn = 0;
    for(j=0;j<64;j++){
      if(j!=i){
	if(Codon[i].aa_index==Codon[j].aa_index){
	  Codon[i].synonym_index[num_syn++]=j;
	}
      }
    }
    Codon[i].num_synonym = num_syn;
  }

  //Define one step synonyms
  for(i=0;i<64;i++){
    num_syn = Codon[i].num_synonym;
    strcpy(codon, Codon[i].codon);
    num_one_step_syn=0;

    //go through each syn
    for(j=0;j<num_syn;j++){
      num_matches=0;

      synon_index = Codon[i].synonym_index[j];
      strcpy(synonym, Codon[synon_index].codon);
      //go through each nt in codon
      for(k=0;k<3;k++){
	if(codon[k]==synonym[k]){
	  num_matches++;
	  codon_match[k]=1;
	}else{
	  codon_match[k]=0;
	}
      }
      if(num_matches==2){
	Codon[i].one_step_synonym_index[num_one_step_syn] = synon_index;
	
	//find mismatch
	k=0;
	while(codon_match[k]==1){
	  k++;
	}
	
	//get nts
	orig_nt =codon[k];
	mut_nt = synonym[k];
	
	//set initial relative mutation rate
	relative_mut_rate = 1;

	switch(mut_nt){
	case 'A':
	case 'a':
	case 'T': 
	case 't':
	  relative_mut_rate*= mu_bias;
	  break;
	  
	default:
	  //Do nothing
	  break;
	  //relative_mut_rate*=(1-at_bias);
	}//end case
	
	//determine if rate is scaled by transition/transversion ratio: gamma_ratio
	if( (orig_nt =='G' && mut_nt =='A') 
	    ||  (orig_nt =='A' && mut_nt =='G') 
	    ||(orig_nt =='C' && mut_nt =='T') 
	    || (orig_nt =='T' && mut_nt =='C')){
	  relative_mut_rate *= gamma_ratio;
	}
	Codon[i].one_step_relative_mutation_rate[num_one_step_syn]= relative_mut_rate;

	num_one_step_syn++; //increase count

      }//end match for one step synonym
 
    }
    Codon[i].num_one_step_synonym = num_one_step_syn;
  }

  //Define factors of (c_i/(c_i+b))/(c_j/(c_j+b)) for each one step synonym
  for(i=0;i<64;i++){
    num_one_step_syn = Codon[i].num_one_step_synonym;
    //go through each one step syn
    for(j=0;j<num_one_step_syn;j++){
      synon_index =  Codon[i].one_step_synonym_index[j];
      Codon[i].elong_pr_ratio[j] = Codon[i].elong_pr/Codon[ synon_index].elong_pr;
    }
  }


  //Define differences of (b/c_i - b/c_j) for each one step synonym
  for(i=0;i<64;i++){
    num_one_step_syn = Codon[i].num_one_step_synonym;
    //go through each one step syn
    for(j=0;j<num_one_step_syn;j++){
      synon_index =  Codon[i].one_step_synonym_index[j];
      Codon[i].delta_b_over_c[j] = B*(1./Codon[i].elong_rate - 1./Codon[synon_index].elong_rate);
    }
  }

}


void Read_Commandline_Args(int argc, char *argv[]){  
  char *dp;
  
  int num_aa, num_codons;
  int i, j, test;
  int outlength;
  int ncheck = 18; //should be equal to size of check[]
  int check[18]={0};// ,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //to add more command line arguments increase this by 0
   //trna, fasta_file, out
  if(print_debug){
    printf("reading arguments\n");
    fflush(stdout);
  }
  for(i=1;i<argc;i++)
    {	if(argv[i][0] == '-')
	{	switch(argv[i][1])
	    {
	    case 'P': //protein production rate scaling factor (phi)
	    case 'p':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nP value not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  check[17]=1;
		  i++;
		  Phi=atof(argv[i]);
		  if(Phi <=0){
		    printf("P value %f is out of range.  Must be greater than 0\n", Phi);
		    wrong();
		  }
		  if(print_debug) printf("Adjusting protein production rate scaling %f\n", Phi);
		  break;
		}



	      
	    case 'T':
	    case 't':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\ntRNA file name not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  //n_aa=Read_Files(argv[++i],1);
		  n_aa=Read_tRNA_File(argv[++i]);
		  strcpy(trna,argv[i]);
		  check[0]=1;
		  break;
		}
	    case 'F':
	    case 'f':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	
		  printf("\nSequence file name not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  //n_seq=Read_Files(argv[++i],2);
		  n_seq=Read_FASTA_File(argv[++i]);
		  strcpy(fasta_file,argv[i]);
		  check[1]=1;
		  if(print_debug){
		    printf("loading file %s\n",fasta_file);
		    fflush(stdout);
		  }
		  break;
		}
	    case 'N':
	    case 'n':
	      {
		switch(argv[i][2])
		  {
		  case 'e':
		  case 'E':
		    if((argv[i][3] != '\0') || (i==argc-1))//check to make sure it's not followed by another letter 
		      //or it's not the last argument
		      {	
			printf("\nIncorrect option or no argument given: %s\nExiting\n", argv[i]);
			wrong();
		      }	
		    else
		      {	
			Ne=atof(argv[++i]);
			break;
		      }
		  default:
		    printf("\nIncorrect argument given: %s\nExiting\n", argv[i]);
		    wrong();
		    break;
		  }
		break;
	      }
	    case 'W':
	    case 'w':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nPrinting file option not specified or Incorrect usage\n");
		  wrong();
		}	
	      else
		{	
		  check[4]=1;
		  i++;
		  pout=atoi(argv[i]);
		  if(print_debug) printf("pout value is %d\n", pout);
		  break;
		}
	    case 'O':
	    case 'o':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nOutput suffix not specified or Incorrect usage\n");
		  wrong();
		}	
	      else
		{	check[5]=1;i++;
		  strcpy(out_prefix,argv[i]);
		  break;
		}
	    case 'A':
	    case 'a':
	      {
		switch(argv[i][2])
		  {
		  case '1':
		    if((argv[i][3] != '\0') || (i==argc-1))
		      {	printf("\nInitiation cost a1 not specified or Incorrect usage\n");
			wrong();
		      }
		    else
		      {	
			check[9]=1;
			i++;
			A1=atof(argv[i]);
			break;
		      }
		  case '2':
		    if((argv[i][3] != '\0') || (i==argc-1))
		      {	printf("\nElongation cost a2 not specified or Incorrect usage\n");
			wrong();
		      }
		    else
		      {	
			check[10]=1;
			i++;
			A2=atof(argv[i]);
			break;
		      }
		  case 'T': //adjust calculations for AT bias. 
		    if((argv[i][3] != '\0') || (i==argc-1))
		      {	printf("\nAT Bias not specified or Incorrect usage. Should be a float\n");
			wrong();
		      }
		    else
		      {	
			//check[11]=1;
			i++;
			at_bias=atof(argv[i]);
			if(at_bias >= 1 || at_bias <= 0){
			  printf("\nAT Bias %f is out of acceptable range. Must be between 0 and 1.\n", at_bias);
			  wrong();

			}
			break;
		      }
		  default: 
		      	printf("\nIncorrect specification of translation costs\n");
			wrong();
		  }
		break;
	      }

	    case 'B':
	    case 'b':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nB value not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  check[11]=1;
		  i++;
		  B=atof(argv[i]);
		  break;
		}

	    case 'D':
	    case 'd':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nD (print out Delta eta) values not specified or Incorrect usage. Should be 0 or 1\n");
		  wrong();
		}
	      else
		{	
		  check[12]=1;
		  i++;
		  print_delta_eta_vals=atoi(argv[i]);
		  break;
		}

	    case 'M':
	    case 'm':
	      //Max time to simulate
	      //If set to a negative value it will be used as a factor times -mu
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nM --max # of time (in generations) to simulate not specified or Incorrect usage.\n");
		  wrong();
		}
	      else
		{	
		  check[13]=1;
		  i++;
		  global_max_time=atof(argv[i]);
		  break;
		}

	    case 'I': //ignore_aa option
	    case 'i':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nI value not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  check[14]=1;
		  i++;
		  ignore_aa=atoi(argv[i]);
		  if(print_debug) printf("ignoring last %d amino acids\n", ignore_aa);
		  break;
		}
	    case 'R': //start at random spot if 1
	    case 'r':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nR value not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  check[15]=1;
		  i++;
		  random_start=atoi(argv[i]);
		  if(print_debug) printf("ignoring last %d amino acids\n", ignore_aa);
		  break;
		}
	    case 'V': //transition to transversion ratio.  
	    case 'v':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nV value not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  check[15]=1;
		  i++;
		  gamma_ratio=atof(argv[i]);
		  if(gamma_ratio <=0){
		    printf("V value %f is out of range.  Must be greater than 0\n", gamma_ratio);
		    wrong();
		  }
		  if(print_debug) printf("Adjusting transition to transversion ratio to %f\n", gamma_ratio);
		  break;
		}
	    case 'Q':
	    case 'q':
	      if((argv[i][2] != '\0') || (i==argc-1))
		{	printf("\nQ value not specified or Incorrect usage\n");
		  wrong();
		}
	      else
		{	
		  check[16]=1;
		  i++;
		  Q=atof(argv[i]);
		  break;
		}
	    }
	}
    }
  // Initializing Defaults, Warnings and Errors
  for(i=0;i<=ncheck;i++)
    { 
      if(print_debug){	printf("processing argument %d \n",i);
	fflush(stdout);
      }

      switch(i){	
      case 0:
	if(check[i]==0)
	  {	
	    dp=&def_tRNA[0];
	    //n_aa=Read_Files(dp,1);
	    n_aa=Read_tRNA_File(dp);
	    strcpy(trna,def_tRNA);
	    break;
	  }
	else if(n_aa<20)
	{
	  printf("\n*Warning*\nt_RNA file contains codons for %d amino acids only\n\n",n_aa);
	  exit(1);
	}
      case 1:
	if(check[i]==0)
	  {	printf("\nSequence file not specified\n");
	    wrong();
	  }
	else if(n_seq<1)
	  {	printf("\nSequence file is empty or not in FASTA format\nCheck README for file formats");
	    wrong();
	  }
      case 2: //carry over from SEMPPR.  Not used here
	if((check[i]==0)||(check[i]<1))
	  {	runs=1000;
	    break;
	  }
      case 3:
	if(check[i]==0)
	  {	MLE=1;
	    break;
	  }
      case 4:
	if(check[i]==0)
	  {	pout=1;
	    break;
	  }
      case 5:
	if(check[i]==0)
	  {	strcpy(out_prefix,"output/out");
	    break;
	  }
      case 6:
	if(check[i]==0)
	  {	analytic=1;
	    break;
	  }
      case 7:
	if(check[i]==0)
	  {	pconf=0;
	    break;
	  }
      case 8:
	if(check[i]==0)
	  {	gconf=0;
	    break;
	  }
      case 9:
	if(check[i]==0)
	  {	A1=4;
	    break;
	  }
      case 10:
	if(check[i]==0)
	  {	A2=4;
	    break;
	  }
      case 11:
	if(check[i]==0)
	  {	B=0.0025;
	    break;
	  }
      case 12:
	if(check[i]==0)
	  {	print_delta_eta_vals=0;
	    break;
	  }
      case 13:
	if(check[i]==0)
	  {	 //global_max_time = 2E10;
	    break;
	  }
      case 14:
	if(check[i]==0)
	  {	ignore_aa = 0;
	    break;
	  }

      }
    }


  //set up out_folder based on out_prefix  
  j = (int)(strlen(out_prefix))-1;
  test = 0;
  //search backwards for '/'
  while(j>0 && test == 0){
    if(out_prefix[j]=='/') {
      test = 1;
    }else{
      j--;
    }
  }

  if(j ==0){ 
    //printf("Printing output to local directory\n");
    strcpy(out_folder, "./");
      }
  else{
    for(i=0;i<j;i++){
      out_folder[i] = out_prefix[i];
    }
    out_folder[j] = '\0';
    if(print_debug){
      printf("using output folder: %s", out_folder);
      fflush(stdout);
    }
  }



  //ensure outfolder exists
  //from man 3p mkdir.  Make with owner and group r/w 
  //mkdir(out_folder, 222); //044 gives d---, 144 gives d-w--w----, 244 gives d-wxrw-r, 344 gives r-x, 224 gives d-wxr-----
  //set runs and MLE to 0 if running in analytic mode

  if(analytic){
    runs=0;
    MLE=0;
  }

//
//  switch(pout)//based on desire output set indicator variables
//    {
//    case -2: //return gmean via stdout
//    case -12: //return id and gmean via stdout
//      calcgmean=1;
//      break;
//    case -11: //return amean via stdout
//    case -1: //return id and amean via stdout
//      calcamean=1;
//      break;
//    case 0: //calc mode-- will be done anyways
//    case -10: //calc mode-- will be done anyways
//      break;
//    default:
//      calcamean = 1;
//      calcgmean = 1;
//      calcvar = 1;
//      break;
//    }      

}



//function to print command line
//taken from K&R p. 115
void Print_Commandline(int argc, char argv[], FILE *outfile){
  while(argc>0)
    fprintf(outfile, "%s%s", *++argv, (argc >1) ? " " : "");
  fprintf(outfile, "\n");
}



void Print_to_Fasta(struct seq_struct *seq,  FILE **outfile){
  int i, j;
  int aa_count;
  int max_chr = 70; //max # of characters/line
  int nt_count;

  aa_count = seq->aa_count;
  fprintf(*outfile, ">%s\tphi=\t%f\n", seq->id, seq->phi_obs);

  nt_count = 0;
  for(i=0;i<=aa_count;i++){//include stop codon
    nt_count+=3;
    if(nt_count<max_chr){
      fprintf(*outfile, "%s", seq->codon_seq[i]);
    }else{
      nt_count-=3;
      j=0;

      while(j<3){
	if(nt_count==max_chr){
	  //start new line
	  fprintf(*outfile, "\n");
	  //reset counter
	  nt_count=0;
	}
	fprintf(*outfile, "%c", seq->codon_seq[i][j]);
	j++;
	nt_count++;
      }
    }
  }
  fprintf(*outfile, "\n");
}
void Old_Print_to_Fasta(struct seq_struct *seq, FILE **outfile){
  int i;
  int aa_count;
  

  

  aa_count = seq->aa_count;
  
  fprintf(*outfile, ">%s\tphi=\t%f", seq->id, seq->phi_obs);
  for(i=0;i<=aa_count;i++){//include stop codon
    if(i%23==0){
      fprintf(*outfile, "\n");
    }
  fprintf(*outfile, "%s", seq->codon_seq[i]);
  }
  fprintf(*outfile, "\n");

}


void Print_Summary_Info(int id_start, int id_stop, int steps)
{
  int i, j, k;

  char filename[150];

  FILE *outfile;
  

  if(id_start <0)
	{
    strcpy(filename, out_prefix);
    strcat(filename,".evol.summary.tsv");
    outfile=fopen(filename,"w+");
	
	fprintf(outfile, "ORF\t");
      fprintf(outfile, "phi\t");
      fprintf(outfile, "Evolve Steps\t");
      fprintf(outfile, "Evolve Time\t");
      fprintf(outfile, "Eta_initial\t");
      fprintf(outfile, "Eta_final\t");
      fprintf(outfile, "Final \\Delta Eta\t");
      fprintf(outfile, "Mean Mutation \\Delta Eta\t");
      fprintf(outfile, "Var Mutation \\Delta Eta\t");
      fprintf(outfile, "# Fix\n");
	id_start=0;
	fclose(outfile);
  }

   strcpy(filename, out_prefix);
    strcat(filename,".evol.summary.tsv");
    outfile=fopen(filename,"a+");

  for(i=id_start; i<=id_stop; i++){
    fprintf(outfile, "%s\t", Seq[i].id);
    fprintf(outfile, "%.4g\t", Seq[i].phi_obs);
    fprintf(outfile, "%d\t", Seq[i].evol_steps);
    fprintf(outfile, "%.0f\t", Seq[i].evol_time);
    fprintf(outfile, "%.3f\t", Seq[i].eta_initial);
    fprintf(outfile, "%.3f\t", Seq[i].eta_obs);
    fprintf(outfile, "%.3f\t", Seq[i].eta_initial - Seq[i].eta_obs);
    fprintf(outfile, "%.4g\t", Seq[i].delta_eta_mean);
    fprintf(outfile, "%.4g\t", Seq[i].delta_eta_var);
    fprintf(outfile, "\n");
  }

  fclose(outfile);
  }


void Print_Output(struct timeval *time_vec, int argc, char *argv[])
{
  int i, j, k;
  int time_counter;
  int start_filename;
  int stop_filename;
  int test;
  FILE *outfile; 
  double read_time, comp_time, print_time, total_time;
  char filename[150];
  char abbr_fasta_file[60];
  char tmpstr[60]; //used for testing end of fasta_filename
  //out_prefix is globally defined




  //Calc read time
  time_counter = 1;
  timersub(&time_vec[time_counter], &time_vec[time_counter-1],&time_diff);
  read_time = (time_diff.tv_sec + time_diff.tv_usec/1000000.);
  
  //Calc computation time
  timersub(&comp_stop,&comp_start,&time_diff);
  comp_time = (time_diff.tv_sec + time_diff.tv_usec/1000000.);


  //print final seq to separate file
  if(pout >=0){ //start printout

    //set prefix
    switch(1){
    case 1:
      strcpy(filename, out_prefix);
      strcat(filename, ".");
      break;
      
    case 2: 
      //put in a separate folder: old approach
      strcpy(filename, out_folder);
      strcat(filename, "/fasta_final/");
      break;

    case 3:  //use a standard name to print output
      strcpy(filename, out_folder);
      i = strlen(filename);
      //check to see if it ends in a /
      //need to offset i by 1 due to C indexing
      if(!(filename[i-1]=='/'))	strcat(filename, "/");
      strcat(filename, "ces.output");

    }

    //decide whether to use input file name for creating output file name
    if(0){
      //process input file name to make output file mimic it in structure
      //get length
      stop_filename = (int)(strlen(fasta_file));
      
      //get last instance of '/'
      i=stop_filename-1;

      while(i>0){
	if('/' ==fasta_file[i]) break;
	i--;
      }
      
      start_filename = ++i;
      for(j=start_filename; j<stop_filename; j++){
	tmpstr[j-start_filename] = fasta_file[j];
      }
      tmpstr[j-start_filename]='\0';
      
      strcat(filename, tmpstr);//append filename
    }

    
    //handle fasta at end
    switch(2){
    case 1:
      //get rid of fasta at end of file name if it exists
      i = (int)(strlen(filename))-6; //move back 6 spots ('fasta/0') 
      for(j=0; j<6; j++){
	tmpstr[j] = filename[i+j];
      }
      tmpstr[j]='\0';
      
      //printf("\n\tPrinting to single file. Fasta file ends with %s\n", tmpstr);
      //test to see if fasta_file name ends with 'fasta', if so truncate
      if(strcmp("fasta", tmpstr)){
	filename[i]='\0';//get rid of '.'
      }
      
      //append to filename
      strcat(filename, ".fasta");
      break;
    
    case 2: //use out prefix and add ces.fasta
      strcat(filename, "ces.fasta");
      break;
    case 3:
      strcat(filename, ".fasta");
      break;
    }

    outfile=fopen(filename, "w+");
    if(outfile==NULL){
      printf("Cannot open file %s. Folder must exist before running program. Exiting.\n", filename);
      exit(1);
    }
  
    fprintf(outfile, "File generated on %s", exec_time);
    fprintf(outfile, "Command Used: \n\t");
    for(i=0; i< argc; i++){
      fprintf(outfile, "%s ", argv[i]);
    }
    fprintf(outfile, "\n");

    for(i=0; i<n_seq; i++){
      Print_to_Fasta(&Seq[i], &outfile);
    }
    fclose(outfile);
  }


  
    gettimeofday(&time_vec[time_counter],NULL);
    timersub(&time_vec[time_counter], &time_vec[time_counter-1],&time_diff);
    print_time = (time_diff.tv_sec + time_diff.tv_usec/1000000.);
    
    timersub(&time_vec[time_counter], &time_vec[0],&time_diff);
    total_time = (time_diff.tv_sec + time_diff.tv_usec/1000000.);

  if(print_debug){
    printf("Calcuated total_time: %f s\n", total_time);
      fflush(stdout);
    }
  
  //  timeval_subtract(&print_time, &print_stop,&comp_stop);
  //timeval_subtract(&total_time, &print_stop,&comp_start);
  //  time(&print_stop);
  //print_time=difftime(print_stop,comp_stop);
  //total_time=difftime(print_stop,start);
  
  //printf("\nCodon Evolution Simulation (CES) Version 0.9\nThe code was executed on %s\nInput file : %s\ntRNA file : %s\nOutput : %s\nNo. of sequences : %d\nNo. of evolutionary steps : %d\nIgnore last AA: %d,\nRuntimes (Read, Computation, Print, Total) : %.4lf\t%.4lf\t%.4lf\t%.4lf seconds\nCommand: %s\n\n\n",exec_time,fasta_file,trna,out_prefix,n_seq, global_max_time,ignore_aa,read_time,comp_time,print_time, total_time, command_line);

  // append log file
  strcpy(filename,out_prefix);
  strcat(filename,".ces.log");
  outfile=fopen(filename,"a+");
  fprintf(outfile,"\nCodon Evolution Simulation (CES) Version 0.9\nThe code was executed on %s\nInput file : %s\ntRNA file : %s\nOutput : %s\nNo. of sequences : %d\nEvolution Run time : %f\nIgnore last AA: %d,\nRuntimes (Read, Computation, Print, Total) : %.4lf\t%.4lf\t%.4lf\t%.4lf seconds\n\n\n",exec_time,fasta_file,trna,out_prefix,n_seq, global_max_time, ignore_aa,read_time,comp_time,print_time, total_time);
  fprintf(outfile, "Command line arguments used to run code: \n\t");
  for(i=0; i< argc; i++){
    fprintf(outfile, "%s ", argv[i]);
  }
  fprintf(outfile, "\n");
  
  fclose(outfile);
  //to get hostname look into hostname.c
  //at  http://www.koders.com/c/fid9A30C5507F7C374C291FEB93CD918B97BEA16C55.aspx?s=md5
  if(print_debug){
    printf("Finished printing to log file\n");
    fflush(stdout);
  }

}



int Evolve_Sequence(struct seq_struct *seq, int aa_count)
{  

  char wt_codon[4];
  char mut_codon[4];
  
  char mut_table_A[3]={'G', 'T', 'C'};
  char mut_table_T[3]={'C', 'A', 'G'};
  char mut_table_G[3]={'A', 'T', 'C'};
  char mut_table_C[3]={'T', 'A', 'G'};
  
  char mut_table_a[3]={'g', 't', 'c'};
  char mut_table_t[3]={'c', 'a', 'g'};
  char mut_table_g[3]={'a', 't', 'c'};
  char mut_table_c[3]={'t', 'a', 'g'};

  char wt_nt;
  char mut_nt;

  unsigned short int D_array[MAX_DDIM][2];//array with codon position and codon index of one step neighbors
  unsigned short int mut_codon_index; 
  unsigned short int wt_codon_index;
  unsigned short int mut_nt_id; //nt that's mutated
  unsigned short int mut_nt_to; //index for mutation of nt to a particular type.
  unsigned short int mut_codon_pos, nt_pos;//position numbers for the codon and nt that's mutated
  unsigned short int wt_aa_index, mut_aa_index;


  int D_dim; //number of dimensions
  int match;
  int fixation;
  int nt_count;
  int i, j, k;
  int evol_steps;
  int synon_mut; 
  int num_fix;
  int num_nonsynon_mut;
  int num_synon_mut;
  
  //indices used for Generate_D_Array
  int Dmax;
  int D_index;
  int codon_index;

  
  //double Ne = Ne;
  double q = Q; //scaling coefficient
  // double mu = MU; //mutation rate per nt
  double phi; //protein production rate
  double qPhi;
  double two_qNePhi;
  double total_mu;//total mutation rate
  double time=0.0;
  double wait_time;
  double expo_param;
  double wt_sigma_vec[MAX_AA];
  double mut_sigma_vec[MAX_AA];
  double factor; //ratio of mut and wt elong_pr, used to rescale sigma_vecs
  double inv_factor; //inverse of above
  double fitness_ratio;
  double fix_pr;
  double eta_obs;
  double evol_time;
  double time_step;
  double wait_parameter;
  double max_time;

  
  //  double sigma_ratio_vec[MAX_AA];
  double b_over_c_vec[MAX_AA];
  double delta_eta_first_term_vec[MAX_AA]; 
  double delta_eta_second_term_vec[MAX_AA]; 

  
  double delta_eta_vec[MAX_DDIM];//array for storing $\Delta \eta_{i,j}$ values for all
  // one step mutants of the resident allele
  double pi_vec[MAX_DDIM];//array for storing values of $\pi(i->j)$ for all one step mutants
  double pi_total; //sum over pi_vec values.
  double mut_pi; //RV representing replacement allele
  double curr_pi; //used to find replacement allele
  double mutation_vec[MAX_DDIM];//array for storing values of mutation rates $\mu_{i->j}$ for all one step mutants
  double pr_vec[MAX_DDIM];//vector of pi * mu values
  double pr_total;//sum over pr_vec values.
  double mut_pr; //RV representing replacement allele
  double curr_pr; //used to find replacement allele


  //double delta_eta_list[MAX_EVOL_STEPS];
  double delta_eta_mean;//mean effect of synon sub for the current wt seq-- def differs from before
  double delta_eta_var;//var effect of synon sub for the current wt seq --definition differs from before
  double delta_eta_sq;//E(\delta \eta^2), used to calc var

  
  double evol_delta_eta_mean;//mean effect of synon sub evolution
  double evol_delta_eta_var;//var effect of synon sub evolution
  double evol_delta_mean;//used to update mean and var
  
  double rand_num;
  double delta_eta; //eta_wt-eta_mut
  double old_sigma_n;
  double old_eta;
  double old_xi;

  double wt_sigma_n;
  double wt_eta;
  double wt_xi;

  double mut_sigma_n;
  double mut_eta;
  double mut_xi;

  double wt_cost;
  double mut_cost;


  double temp;
  struct tm * timeinfo;
  struct timeval curr_time;

  const gsl_rng_type * rng_type;
  gsl_rng * rn;  

  char delta_eta_filename[150];
  FILE *delta_eta_outfile;


  char fasta_filename[150];
  FILE *fasta_outfile;


  max_time = seq->max_time;
  //open files for writing output





  if(print_evol_fasta){
    strcpy(fasta_filename, out_folder);
    strcat(fasta_filename, "/fasta_evol");

    //make dir to ensure it exists
    //mkdir(fasta_filename, 664);//S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    strcat(fasta_filename, "/");
    strcat(fasta_filename, seq->id);
    strcat(fasta_filename,".evol.fasta");
    fasta_outfile = fopen(fasta_filename, "w+");
  }


  if(print_delta_eta_vals){
    strcpy(delta_eta_filename, out_folder);
    strcat(delta_eta_filename, "/delta_eta_lists");

    //make dir to ensure it exists
    //mkdir(delta_eta_filename, 664);//S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


    strcat(delta_eta_filename, "/");
    strcat(delta_eta_filename, seq->id);
    strcat(delta_eta_filename,".delta_eta_list.tsv");
    delta_eta_outfile = fopen(delta_eta_filename, "w+");
    if(delta_eta_outfile==NULL){
      printf("Cannot open file %s. Folder must exist before running program. Exiting.\n", delta_eta_filename);
      exit(1);
    }

  }

  //Set up GSL RNG
  gsl_rng_env_setup();
  

  
  //Create a generator chosen by the 
  //environment variable GSL_RNG_TYPE
  //use default gsl for generating uniform rn 
  //from which other rn functions are derived 

  rng_type = gsl_rng_default;
  rn = gsl_rng_alloc (rng_type);

  //seed rng using clock
  gettimeofday(&curr_time,NULL);
  if(print_debug){
    printf("Seeding RNG with clock. Current time is %d %d\n", (int)(curr_time.tv_sec), (int)(curr_time.tv_usec));
    fflush(stdout);
  }

  gsl_rng_set(rn, (curr_time.tv_sec + curr_time.tv_usec));// (const gsl_rng * r, unsigned long int s)   

  //end RNG set up


  phi = seq->phi_obs;
  qPhi = Q*phi;
  two_qNePhi = qPhi* 2* Ne;
  

  //references C++
  //two ways of passing the seq.sigma_vec
  //&seq->sigma_vec
  //&(*seq).sigma_vec
  //However, these seem to pass a pointer to an array as opposed to a pointer
  //to the start of the array, so the syntax that seems to work for my purposes is
  // &(seq->codon_index[0]) gives the exact address I want.

  //seq->value would work for a scalar

  //  sigma_vec = &seq->sigma_vec;

  Convert_Sigma_Vec_to_Sigma_Ratio_Vec(&(seq->sigma_vec[0]), &(seq->sigma_ratio_vec[0]), aa_count);

  //calculate b_over_c_vec
  Calc_B_over_C_Vec(&(seq->codon_index[0]), b_over_c_vec, aa_count);

  //end Generate_D_Array
  

  evol_steps = 0;
  /////////EVOLUTION OF SEQUENCE USED TO OCCUR HERE
      for(evol_time=0; evol_time < max_time; evol_time+=time_step){
    // simulation begins.  Ideally could start further down 
    // and just update below terms in just a few places
	
	evol_steps++;
	


    //determining Dimensionality/# one step neighbors of the system.
    D_dim = Calc_Dimensionality(&(seq->codon_index[0]), aa_count);
    

    //Create D_index of position and index values 
    //Generate_D_Array(seq->codon_index, &D_array, aa_count);
    //  imax = (unsigned short int)(aa_count);
    D_index=0;
    for(i = 0;i<aa_count;i++){
      codon_index =seq->codon_index[i];
      Dmax = Codon[codon_index].num_one_step_synonym;
      
      //go through each synonym
      for(j=0;j<Dmax;j++){
	//Get mutation rate
	mutation_vec[D_index] = Codon[codon_index].one_step_relative_mutation_rate[j];
			 
        (D_array[D_index][0])= (unsigned short int)(i);
	(D_array[D_index++][1])= (Codon[codon_index].one_step_synonym_index[j]);


      }
      
    }
    
    if(print_debug && 1==1){
      printf("Step: %d: Time of Replacement = %.0f\n\tCalculating Delta_Eta First and Second Terms\n", evol_steps, evol_time+time_step);
      fflush(stdout);
    }

    Calc_Delta_Eta_First_and_Second_Term_Vecs(&(seq->sigma_ratio_vec[0]), b_over_c_vec, delta_eta_first_term_vec, delta_eta_second_term_vec, aa_count);
    
    if(print_debug && 1==1){
      printf("\tCalculating Delta Eta Vec\n");
      fflush(stdout);
    }

    Calc_Delta_Eta_Vec(&(seq->codon_index[0]), b_over_c_vec, delta_eta_first_term_vec, delta_eta_second_term_vec, delta_eta_vec, aa_count, D_dim);
    
    
    if(print_debug && 1==1){
      printf("\tCalculating Pi Vec\n");
      fflush(stdout);
    }

    //Calculate pr of allele replacement given a mutant arises
    Calc_Pi_Vec(delta_eta_vec, pi_vec, &pi_total, phi, qPhi, two_qNePhi, aa_count, D_dim);

    //test to see if there are mutation effects
    if((mu_bias==0.5) && gamma_ratio ==1.){
      //flat mutation so pi values = pr values
      //Assign pr_vec address to pi_vec address
      //This should work so long as the dimensions of these
      //vectors are the same
      //      &(pr_vec[0]) = &(pi_vec[0]);//way of copying pointers according to Kernighan and Ritchie p. 98
      //The above doesn't work because pr_vec represents a specific space in memory and it its address can't
      //be reassigned.  It can be copied, however.

      memcpy(pr_vec, pi_vec, D_dim*sizeof(double));

      pr_total = pi_total;

    }else{
      //Calculate true pr of replacement, weighting by mutation rates
      Calc_Replacement_Pr_Vec(mutation_vec, pi_vec, pr_vec, &pr_total, D_dim);
    }


    //Calculate expected time to mutation to allele that will replace resident
    //We are assuming that replacement happens quickly.
    wait_parameter = 1/(Ne * mu * pr_total);// This is the 'failure' or leaving rate of the resident allele.

    //ideally we would expect the wait parameter to be >>1, but it likely doesn't matter 
    //if(wait_parameter < 100)

    //The time to replacement should follow a geometric dist which we approximate 
    // with an exponential distribution.
    //Note the parameter for GSL's RNG argument is the expected wait time, so it is 1/replacement rate
    time_step = gsl_ran_exponential(rn, wait_parameter);


    mut_pr = gsl_rng_uniform(rn) * pr_total;
    if(print_delta_eta_vals){
      if(evol_time==0){
	//print header
	fprintf(delta_eta_outfile, "evol_step\tevol_time\tSum_Fix_Pr\teta_wt");
	for(i=0; i<D_dim; i++){
	  fprintf(delta_eta_outfile, "\tdelta_eta_%d", (i+1));
	}
	fprintf(delta_eta_outfile, "\n");
      }
      fprintf(delta_eta_outfile, "%d\t%.6e\t%.6e\t%.6e", evol_steps, evol_time+time_step, pr_total, wt_eta);
      for(i=0; i<D_dim; i++){
	fprintf(delta_eta_outfile, "\t%.6e", delta_eta_vec[i]);
      }
      fprintf(delta_eta_outfile, "\n");
      
    }


    //evaluate moments of delta_eta_vec
    delta_eta_mean = 0.0;
    delta_eta_sq = 0.0;
    for(i=0;i<D_dim;i++){
      delta_eta_mean +=delta_eta_vec[i];
      delta_eta_sq+=pow(delta_eta_vec[i], 2);
    }

    delta_eta_mean/=D_dim;
    delta_eta_sq/=D_dim;
    
    delta_eta_var = delta_eta_sq - pow(delta_eta_mean, 2);

    seq->delta_eta_mean = delta_eta_mean;
    seq->delta_eta_var = delta_eta_var;

    //it seems like for(i=1, i< max_time, i++){ loop could start here if I was to 
    //think about it more.
    // The problem is D_Array set up kinda screwy so that one may need to shift
    // everything up or down depending on the # of neighbors of the mutant
    // Further, lots of stuff needs to be rescaled
    // 
    if(print_debug && 1==1){
      printf("\tChoosing mutant...");
      fflush(stdout);
    }
    
    //choose a RV to determine which allele replaces the resident
    mut_pr = gsl_rng_uniform(rn) * pr_total;

    //Note: pr_total is not scaled by the constants Ne or mu)
    if(print_debug && 1==1){
      printf("%g from 0 to %g \n", mut_pr, pr_total);
      fflush(stdout);
    }
  
    //Find out which allele it is by searching forwards. 
    //This is likely more efficient than starting at the end b/c there is weaker selection at 
    // the start of the sequence
    curr_pr = 0.0;
    i = 0;
    while(curr_pr < mut_pr){
      curr_pr+=pr_vec[i++];
    }
    
    i--;
    D_index = i;
    
    
    mut_codon_pos = D_array[D_index][0];
    mut_codon_index = D_array[D_index][1];

    //get wt codon index
    wt_codon_index = seq->codon_index[mut_codon_pos];

    //get codons
    strcpy(wt_codon,Codon[wt_codon_index].codon);
    strcpy(mut_codon,Codon[mut_codon_index].codon);


    delta_eta = delta_eta_vec[i];//should be \eta_i - \eta_j
    
    //values to check to make sure delta_eta is working right
    if(1==0){
      wt_xi = Calc_Xi(&(seq->sigma_vec[0]), aa_count);
      wt_sigma_n = seq->sigma_vec[aa_count-1];
      wt_eta = Calc_Eta(wt_sigma_n, wt_xi, aa_count);
      seq->eta_initial = wt_eta;
    }

    //update everything//

    //update evol_delta_mean/var vals
    evol_delta_mean = delta_eta -delta_eta_mean;
    evol_delta_eta_mean +=evol_delta_mean/(evol_time+time_step);
    evol_delta_eta_var += evol_delta_mean * (delta_eta - evol_delta_eta_mean);
    //update values for printing
    seq->evol_delta_eta_mean = evol_delta_eta_mean;
    seq->evol_delta_eta_var = evol_delta_eta_var;

    //update codon 
    strcpy(seq->codon_seq[mut_codon_pos],Codon[mut_codon_index].codon);

    //update codon index
    seq->codon_index[mut_codon_pos] = mut_codon_index;

    
    //update Sigma_Ratio_Vec by dividing the values up to mut_codon_pos by elong_pr_i/elong_pr_j
    //values after this point don't change since they involve the factor in both the
    //denominator and numerator

    factor = Codon[mut_codon_index].elong_pr/Codon[wt_codon_index].elong_pr;
    inv_factor = Codon[wt_codon_index].elong_pr/Codon[mut_codon_index].elong_pr;
  
    for(i=1;i<=mut_codon_pos;i++){
      seq->sigma_ratio_vec[i]*=inv_factor;
    }

    //update sigma_vec for values at and above mut_codon_pos;
    for(i=mut_codon_pos; i< aa_count; i++){

      seq->sigma_vec[i]*=factor;
    }

    //update b_over_c_vec
    b_over_c_vec[mut_codon_pos] = B/Codon[mut_codon_index].elong_rate;

    //update D_array
    //    D_array[D_index][1]=wt_codon_index;


    //update eta_obs by subtracting difference: \eta_wt - (\eta_wt - \eta_mut) = \eta_mut
    seq->eta_obs-=delta_eta;


    if(1==0){
      //check to make sure delta_eta is working right
      mut_xi = Calc_Xi(&(seq->sigma_vec[0]), aa_count);
      mut_sigma_n = seq->sigma_vec[aa_count-1];
      mut_eta = Calc_Eta(mut_sigma_n, mut_xi, aa_count);
      
      if(abs((wt_eta - mut_eta)- delta_eta)/delta_eta > 1.0E-6){
	
	printf("ERROR: wt_eta-mut_eta = %f != %f delta_eta\n", (wt_eta - mut_eta), delta_eta);
      }
    }

    //    eta_obs +=delta_eta;

    
    if(print_debug || 1==0){
      temp = (evol_time+time_step);
      if(evol_steps ==1) printf("\tevol_time\tmut_codon_pos\twt_codon\tmut_codon\tdelta_eta\tdelta_eta_mean\tdelta_eta_var\n");
      printf("\t%g\t%d\t%s\t%s\t%g\t%g\t%g\n", temp, mut_codon_pos, wt_codon, mut_codon, delta_eta, seq->delta_eta_mean, seq->delta_eta_var);
      
      
    }
  } //end sequence evolution

      

  if(print_evol_fasta){
    fclose(fasta_outfile);
  }

  if(print_delta_eta_vals){
    fclose(delta_eta_outfile);
  }

  if(pout <= 0){ //print to individual fasta files for each gene
  //print final seq to separate file
  strcpy(fasta_filename,out_folder);
  strcat(fasta_filename, "/fasta_final/");

  //make dir to ensure it exists
  //mkdir(delta_eta_filename, 664);//S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
  strcat(delta_eta_filename, "/");

  strcat(fasta_filename, seq->id);
  strcat(fasta_filename, ".final.fasta");

  fasta_outfile=fopen(fasta_filename, "w+");
  if(fasta_outfile==NULL){
    printf("Cannot open file %s. Folder must exist before running program. Exiting.\n", fasta_filename);
    exit(1);
  }
  
 
    Print_to_Fasta(seq, &fasta_outfile);
    fclose(fasta_outfile);
  }

  //calc mean
//
//  delta_eta_mean = 0;
//  for(i=0; i<num_fix;i++){
//    delta_eta_mean+=delta_eta_list[i];
//  }
//  delta_eta_mean /=num_synon_mut;
//
//  //calc var
//  delta_eta_var = 0;
//  for(i=0; i<num_fix;i++){
//    delta_eta_var+=pow(delta_eta_list[i], 2);
//  }
//  delta_eta_var /=(num_synon_mut-1);

  //update Seq[]
  seq->evol_steps =evol_steps; 
  seq->evol_time =evol_time; 
  seq->delta_eta_mean = delta_eta_mean;
  seq->delta_eta_var = delta_eta_var;
  

  // printf("%s mutations with %d synon and  %d fixations. \nRatios %f and %f\n\n", max_time,num_synon_mut, num_fix, (double)(num_fix)/max_time, (double)(num_fix)/num_synon_mut);
  gsl_rng_free (rn);

  return 1;
}






// Main Function
int main(int argc, char *argv[])
{	
  int i, j, k, D_dim;
  int time_counter, aa_count;
  double max_time;
  time_t start;
  double time_spent[5]; //time spent in each step
  struct tm * timeinfo;
  struct timeval time_val, time_diff;//, comp_start, comp_stop, print_stop, time_diff;
  struct timeval time_vec[5]; //0 prog_start, 1 read_stop, 2 comp_stop, 3 print_stop
  double prcnt=0.9999;
  char tmp_str[100];
  time_t time_raw_format;
    
  //get run start time
  time( &time_raw_format );
  sprintf(exec_time, "%s", ctime(&time_raw_format));



  // Initializing the random number generator
  srand (time(NULL));//should this be NULL? mikeg
	
  time_counter=0;
  time(&start);
  timeinfo=localtime(&start);

  gettimeofday(&time_vec[time_counter++],NULL);

  //read in any arguments passed on the command line
  Read_Commandline_Args(argc, argv);


  //save command line info to a string
  //based on code from K&R p. 115

  command_line[0]='\0';
  for(i=0; i<argc; i++){
    sprintf(tmp_str, "%s%s", argv[i], (i<argc-1) ? " " : "");
    strcat(command_line, tmp_str);
  }

  //adjust aa_counts
  if(ignore_aa > 0){
    for(i=0; i<n_seq; i++){
      Seq[i].aa_count-=ignore_aa;
    }
  }

  gettimeofday(&time_vec[time_counter],NULL);
  time_counter++;


  Process_AA_Information();
  //initialize codon structures


  Generate_Codon_Structures();
  

  // Main code begns
  Print_Summary_Info(-1, -1, 0);


  if(global_max_time <0){
    max_time = -global_max_time/mu * 1.0;
    
  }else{
    max_time =  global_max_time *1.0;
  }
  
#pragma omp parallel for private(i, D_dim,aa_count) shared(Seq,AA,Codon, max_time)
  //"pragma omp parallel" is the command which will parallelize the next full command (here a for loop)
  // private is defined for individual processes inaccessible by others
  // shared is for variables accessible to all processes
  
  for(i=0;i<n_seq;i++)
    {//parallelized for loop
      aa_count = Seq[i].aa_count;

      //choose outcome
      if(print_debug || 1==0){
	printf("%d:%s\tProcessing\n", i, Seq[i].id);
	fflush(stdout);
      }


      if(print_debug){
	printf("%d:%s\tConverting Codon strings to index values\n", i, Seq[i].id);
	fflush(stdout);
      }
      Convert_Codon_Seq_to_Codon_Index(Seq[i].codon_seq, Seq[i].codon_index, aa_count);

      if(print_debug){
	printf("Calcuating Dimensionality in main() as %d\n",D_dim);
	fflush(stdout);
      }

      //determining Dimensionality/# one step neighbors of the system.
      D_dim = Calc_Dimensionality(Seq[i].codon_index, aa_count);

      if(print_debug){
	printf("%d:%s\tConverting Codon index values to amino acid index values\n", i, Seq[i].id);
	fflush(stdout);
      }
      Convert_Codon_Index_to_AA_Index(Seq[i].codon_index, Seq[i].aa_index,aa_count);

      if(random_start){
	//replace codn sequence and index  with a randomly generated one
	Generate_Random_Codon_Seq_and_Index(Seq[i].codon_seq, Seq[i].codon_index, Seq[i].aa_index, aa_count);
      }
      
      if(print_debug){
	printf("\t\t\tCalculating Sigma Vec\n", i, Seq[i].id);
	fflush(stdout);
      }
      Convert_Codon_Index_to_Sigma_Vec(Seq[i].codon_index, Seq[i].sigma_vec, aa_count);

      Seq[i].sigma_obs = Seq[i].sigma_vec[aa_count-1];

      
      if(print_debug){
	printf("\t\t\tCalculating Xi, and Eta\n");
	fflush(stdout);
      }
      

      Seq[i].xi_obs = Calc_Xi(Seq[i].sigma_vec, aa_count);

      Seq[i].eta_initial =Seq[i].eta_obs=Calc_Eta(Seq[i].sigma_obs, Seq[i].xi_obs, aa_count);
      
      Seq[i].max_time = max_time;

      Evolve_Sequence(&Seq[i], aa_count);

      Print_Summary_Info(i, i, max_time);

      //      Old_Print_to_Fasta(i);

      

    }//end parallel part of code

  if(print_debug){
    printf("Finished calculations... printing output..\n");
    fflush(stdout);
  }


  gettimeofday(&time_vec[time_counter],NULL);


  Print_Output(time_vec, argc, argv);
	
  
}

