#include "read_file.h" //used to read in fasta and trna files
#include <iostream>
int main(int argc, char *argv[]){
  int aa_count;
  
    //adjust aa_counts
  if(ignore_aa > 0){
    for(int i=0; i<n_seq; i++){
      Seq[i].aa_count-=ignore_aa;
    }
  }
  Read_Commandline_Args(argc, argv);
  for (int i=0; i<n_seq; i++){
    Convert_Codon_Seq_to_Codon_Index(Seq[i].codon_seq, Seq[i].codon_index, aa_count);
    Convert_Codon_Index_to_AA_Index(Seq[i].codon_index, Seq[i].aa_index,aa_count);
    for (int j=0; j<MAX_AA; j++){
      for (int k=0; k<4; k++){
	//	std::cout<<Seq[i].codon_seq[j][k];
	//	std::cout<<Seq[i].aa_index[j];
	std::cout<<Seq[i].codon_index[0];
      }
    }
    std::cout<<std::endl<<std::endl;
  }
    /*

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
    */
}
    
