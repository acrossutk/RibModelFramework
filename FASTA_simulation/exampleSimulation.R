library(AnaCoDa)
library(seqinr)
## FASTA file with sequences you want to simulate
genome.nosp <- initializeGenomeObject("Data/FONSE/nse2000.fasta")

##File containing phi values, not on log scale
## pass in as vector to initial.expression.value in initializeParameterObject, should be in same order as genes in genome
phi <- read.table("Data/FONSE/nse2000.phi.csv",sep=",",header=T)
#Not sure about these
sel  <- "Data/FONSE/selection2ref.csv"
mut <-  "Data/FONSE/nse2000.logmu.csv"

#sel  <- "../Results/mp_sp_pseudo_sanity_check/chain/run_15/Parameter_est/sp_main_Selection"
#mut <-  "../Results/mp_sp_pseudo_sanity_check/chain/run_15/Parameter_est/sp_main_Mutation"
parameter.nosp <- initializeParameterObject(genome.nosp,sphi=c(0.01),num.mixtures = 1,gene.assignment = rep(1,length(genome.nosp)),mixture.definition = "allUnique",initial.expression.values =10^phi[,2],split.serine = TRUE)

## Pass in files amd mixture category, 1 == first mixture
parameter.nosp$initSelectionCategories(sel,1)
parameter.nosp$initMutationCategories(mut,1)
model.nosp <- initializeModelObject(parameter.nosp,model="FONSE",with.phi=FALSE)
model.nosp$simulateGenome(genome.nosp)
genome.nosp$writeFasta(paste0("Data/Simulated/sim1",".fasta"),simulated=T)



#genome.nosp$writeFasta(paste0("../Data/Genomes/Simulated/Simulated_nosp/nosp_sim_sp_dEta",i,".fasta"),simulated=T)

