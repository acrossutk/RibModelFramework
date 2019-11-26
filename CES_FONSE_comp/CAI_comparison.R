
library(AnaCoDa)
library(seqinr)
#will get the CAI for CES and FONSE and compare them using cor() function (correlation coefficient) for all genes in a genome
CAI_corr_coeff <- function(FASTA_ref.nosp,FASTA_CES.nosp, FASTA_FONSE.nosp) {
   CAI_CES <- getCAI(referenceGenome = FASTA_ref.nosp, testGenome = FASTA_CES.nosp)
  CAI_FONSE <- getCAI(referenceGenome = FASTA_ref.nosp,testGenome = FASTA_FONSE.nosp)
return(cor(CAI_CES,CAI_FONSE))
}

a1_data <- data.frame(a1f=double(),correlationf=double())
for (a1 in c(1,2,3,4,6,8,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400)){
filename_CES<-paste("CES_FASTA/S.cerevisiae.S288_a=",a1,".beyer.phi.fasta",sep="")
filename_FONSE<-paste("CES_FASTA/S.cerevisiae.S288_a=",a1,".beyer.phi.fasta",sep="")

FASTA_ref.nosp <- initializeGenomeObject("REF_FASTA/S.cerevisiae.S288.beyer.phi.fasta")
FASTA_CES.nosp <- initializeGenomeObject(filename_CES)
FASTA_FONSE.nosp <- initializeGenomeObject(filename_FONSE)

correlation<-CAI_corr_coeff(FASTA_ref =FASTA_ref.nosp, FASTA_CES = FASTA_CES.nosp,FASTA_FONSE = FASTA_FONSE.nosp)
new_row<-data.frame(a1f=a1,correlationf=correlation)
a1_data<-rbind(a1_data,new_row)
}

