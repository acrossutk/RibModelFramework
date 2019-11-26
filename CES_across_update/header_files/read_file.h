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
#include "global_variables.h" //list of globally defined variables
/*
extern struct Amino_acid;
extern struct Amino_acid AA[22];
extern struct seq_history_struct;
extern struct seq_struct Seq[MAX_LOCI];
*/


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
      Seq[i].phi_obs = strtod(char_ptr, NULL);//first argument is ptr to str location of double
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


void Read_Commandline_Args(int argc, char *argv[]){  
  char *dp;
  
  int num_aa, num_codons;
  int i, j, test;
  int outlength;
  int ncheck = 17; //should be equal to size of check[]
  int check[17]={0}; //,0,0,0,0,0,0,0,0,0,0,0,0,0}; //to add more command line arguments increase this by 0
   //trna, fasta_file, out
  if(print_debug){
    printf("reading arguments\n");
    fflush(stdout);
  }
  for(i=1;i<argc;i++)
    {	if(argv[i][0] == '-')
	{	switch(argv[i][1])
	    {	
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
void Convert_Codon_Index_to_AA_Index(unsigned short int *codon_index_vec, unsigned short int *aa_index_vec, int aa_count)
{
  int i, j, k, match;
  for(i=0; i< aa_count; i++){
    aa_index_vec[i]=Codon[codon_index_vec[i]].aa_index;
  }

}
