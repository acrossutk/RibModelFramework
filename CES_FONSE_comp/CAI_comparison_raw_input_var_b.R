
library(AnaCoDa)
library(seqinr)
#will get the CAI for CES and FONSE and compare them using cor() function (correlation coefficient) for all genes in a genome
CAI_corr_coeff <- function(FASTA_ref.nosp,FASTA_CES.nosp, FASTA_FONSE.nosp) {
  CAI_CES <- getCAI(referenceGenome = FASTA_ref.nosp, testGenome = FASTA_CES.nosp)
  print(CAI_CES)
  CAI_FONSE <- getCAI(referenceGenome = FASTA_ref.nosp,testGenome = FASTA_FONSE.nosp)
return(cor(CAI_CES,CAI_FONSE))
}
Add_phi_values <- #adds phi values to fasta file for CES


b=0.00025
runtime=-50
population=1
filename_CES<-paste("CES_output/vary_b/b=",b,".ces.fasta",sep="")
#CES
runtime_in<-paste('-M',runtime,sep=' ')
b_in<-paste('-B',b,sep=' ')
population_in<-paste('-Ne',population,sep=' ')
FASTA_in<-paste('-F', 'REF/S.cerevisiae.S288c.fasta',sep=" ")
tRNA_in<-paste('-T', 'REF/S.cerevisiae.PNAS.elongation.rates.tsv',sep=" ")
mutation_in=paste('-P', 'REF/S.cerevisiae.PNAS.mutation.rates.tsv',sep=" ")
filename_CES_in<-paste('-O',"CES_FASTA/vary_b/b=",b,".ces.fasta",sep="")

raw_input<-read.table(file = 'REF/S.cerevisiae.S288c.fasta')

CES_system_input<-paste('./../CES_1.0/source/CES',runtime_in,b_in,population_in,FASTA_in,tRNA_in,mutation_in,filename_CES_in,sep=" ")
system(CES_system_input)

#ANaCoDa
phi <- 'REF/S.cerevisiae.beyer.phi.rescaled.tsv'
sel  <- "Data/FONSE/Scereviciae.sel.csv"
mut <-  "Data/FONSE/Scereviciae.mut.csv"
parameter.nosp <- initializeParameterObject(genome.nosp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.nosp)),mixture.definition = "allUnique",split.serine = TRUE, model="FONSE")



#b_data <- data.frame(bf=double(),correlationf=double())
#for (b in c(0.0000025, 0.000025, 0.00025, 0.0025, 0.025)){
#filename_CES<-paste("CES_FASTA/vary_b/b=",b,".ces.fasta",sep="")
#filename_FONSE<-paste("FONSE_FASTA/vary_b/S.cerevisiae.S288_b=",b,".beyer.phi.fasta",sep="")

#FASTA_ref.nosp <- initializeGenomeObject("REF_FASTA/S.cerevisiae.S288.beyer.phi.fasta")
#FASTA_CES.nosp <- initializeGenomeObject(filename_CES)
#FASTA_FONSE.nosp <- initializeGenomeObject(filename_FONSE)

#correlation<-CAI_corr_coeff(FASTA_ref =FASTA_ref.nosp, FASTA_CES = FASTA_CES.nosp,FASTA_FONSE = FASTA_FONSE.nosp)
#new_row<-data.frame(bf=b,correlationf=correlation)
#b_data<-rbind(b_data,new_row)
#}

#b_data <- data.frame(bf=double(),correlationf=double())
#b=0.00025
#  filename_CES<-paste("CES_FASTA/vary_b/b=",b,".ces.fasta",sep="")
#  filename_FONSE<-paste("FONSE_FASTA/vary_b/S.cerevisiae.S288_b=",b,".beyer.phi.fasta",sep="")
  
#  FASTA_ref.nosp <- initializeGenomeObject("REF_FASTA/S.cerevisiae.S288.beyer.phi.fasta")
#  FASTA_CES.nosp <- initializeGenomeObject(filename_CES)
#  FASTA_FONSE.nosp <- initializeGenomeObject(filename_FONSE)
  
#  correlation<-CAI_corr_coeff(FASTA_ref =FASTA_ref.nosp, FASTA_CES = FASTA_CES.nosp,FASTA_FONSE = FASTA_FONSE.nosp)
#  new_row<-data.frame(bf=b,correlationf=correlation)
#  b_data<-rbind(b_data,new_row)


#  CAI_CES <- getCAI(referenceGenome = FASTA_ref.nosp, testGenome = FASTA_CES.nosp)
#  print(CAI_CES)
#  write.table(CAI_CES, "CES.txt")
#  CAI_FONSE <- getCAI(referenceGenome = FASTA_ref.nosp,testGenome = FASTA_FONSE.nosp)
#  print(CAI_FONSE)
#  write.table(CAI_FONSE,"FONSE.txt")
#  cor(CAI_CES,CAI_FONSE)
  