fasta_input<-read.table(file = 'REF/S.cerevisiae.S288c.fasta')
phi_input<-read.table(file='REF/S.cerevisiae.beyer.phi.rescaled.tsv')

phi_names=row.names(phi_input)
phi_names<-paste('>',phi_names,sep="")
#fasta_name=">YAL002W"


for (fasta_row in 1:nrow(fasta_input)){
  
 # if (grepl(">",fasta_input[fasta_row,1],TRUE)==TRUE){
#    print(fasta_input[fasta_row,1])
#  }
  
  if (fasta_row+1 % 2==0){
    
    print(fasta_input[fasta_row,1])
  }
 #   print(fasta_input[fasta_row]) 
  
  
#  for (phi_row in 1:nrow(phi_input)){
#    if (phi_names[phi_row] == fasta_name){
#      print(phi_names[row])
#    }
#  }
}
