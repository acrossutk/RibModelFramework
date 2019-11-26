options(scipen=999)
library(AnaCoDa)
library(seqinr)
## FASTA file with sequences you want to simulate
genome.nosp <- initializeGenomeObject("S.cerevisiae.S288.beyer.phi.fasta")

##File containing phi values, not on log scale
## pass in as vector to initial.expression.value in initializeParameterObject, should be in same order as genes in genome
phi <- read.table("Data/FONSE/Scereviciae.phi.csv",sep=",",header=T)
sel  <- "Data/FONSE/Scereviciae.sel.csv"
mut <-  "Data/FONSE/Scereviciae.mut.csv"

for (b in c(0.0000025, 0.000025, 0.00025, 0.0025, 0.025)){
  print(b)
#parameter.nosp <- initializeParameterObject(genome.nosp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.nosp)),mixture.definition = "allUnique",initial.expression.values =10^phi[,2],split.serine = TRUE)
parameter.nosp <- initializeParameterObject(genome.nosp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.nosp)),mixture.definition = "allUnique",split.serine = TRUE, model="FONSE")
## Pass in files amd mixture category, 1 == first mixture
parameter.nosp$initSelectionCategories(sel,1)
parameter.nosp$initMutationCategories(mut,1)
model.nosp <- initializeModelObject(parameter.nosp,model="FONSE",with.phi=TRUE)
model.nosp$simulateGenome(genome.nosp)
filename<-paste("S.cerevisiae.S288_b=",b,".beyer.phi.fasta",sep="")
genome.nosp$writeFasta(paste0("simulated_data/vary_b/",filename),simulated=T)
}


#genome.nosp$writeFasta(paste0("../Data/Genomes/Simulated/Simulated_nosp/nosp_sim_sp_dEta",i,".fasta"),simulated=T)

