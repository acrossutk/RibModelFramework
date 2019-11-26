library(ribModel)
rm(list=ls())
#read genome

#Initializes the files used for the following parameters
genome.file <- "simulatedFONSE_b-0.00515.fasta"
phi.file <- "simulatedFONSE.phi.csv"
mut.file <- "S.cer.mut.ref.csv"
sel.file <- "deltaomega_b-0.00515.csv"

run.name <- "simulatedFONSE_b-0.001_scaled10000"

from.good.values <- TRUE
genome <- initializeGenomeObject(genome.file)

#initialize sphi and the number of mixtures, but with FONSE, the 
#number of mixtures is always 1. Then initialize the parameter
sphi_init <- 1.2
numMixtures <- 1
mixDef <- "allUnique"
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457), rep(1, 3903))
#geneAssignment <- c(rep(1,448), rep(1,513), rep(2,457))
#geneAssignment <- c(rep(1,448), rep(2,457))
#geneAssignment <- c(rep(1,500), rep(2,500))
geneAssignment <- rep(1,length(genome))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment, model="FONSE", split.serine = TRUE,
                                       mixture.definition = mixDef)

# initialize MCMC object. The MCMC loop will go for samples * thinning iterations.
samples <- 1000
thinning <- 100
adaptiveWidth <- 10
mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth, 
                             est.expression=TRUE, est.csp=TRUE, est.hyper=TRUE)
# get model object
model <- initializeModelObject(parameter, "FONSE")

#TODO: ask Jeremy about the purpose of this boolean check
if(from.good.values) {
  parameter$initMutationCategories(c(mut.file), 1)
  parameter$initSelectionCategories(c(sel.file), 1)
  phi.values <- parameter$readPhiValues(phi.file)
  
  parameter$initializeSynthesisRateByList(phi.values)
}

#Set the settings for the restart file, i.e. the names of them and how often to write them.
setRestartSettings(mcmc, "restartFile.rst", 100, TRUE)

#run mcmc on genome with parameter using model
system.time(
  runMCMC(mcmc, genome, model, 16)
)

full.name <- paste(run.name, samples*thinning, sep="_")

#Write the .Rdat files for both the parameter object and the MCMC object.
writeParameterObject(parameter, file=paste(full.name, "ParamObject.Rdat", sep=""))
writeMCMCObject(mcmc, file=paste(full.name, "MCMCObject.Rdat", sep=""))

#This line makes it so that the plots of the traces go to a pdf instead of being shown in RStudio
pdf(paste(full.name, "MCMC.pdf", sep=""))

#Plots the trace of the log likelihood
plot(mcmc)
loglik.trace <- mcmc$getLogLikelihoodTrace()
acf(loglik.trace)

# plots different aspects of trace. Mixture probability can be ignored since the mixture is always one.
# Also need to ask Jeremy why SPhi doesn't plot anything.
trace <- parameter$getTraceObject()
plot(trace, what = "MixtureProbability")
plot(trace, what = "SPhi")
plot(trace, what = "ExpectedPhi")

mixtureAssignment <- unlist(lapply(1:length(genome),  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples*0.1, geneIndex)}))
expressionValues <- unlist(lapply(1:length(genome), function(geneIndex){
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixtureAssignment[geneIndex])
  parameter$getSynthesisRatePosteriorMeanByMixtureElementForGene(samples*0.1, geneIndex, expressionCategory)
}))
expressionValues <- log10(expressionValues)
obs.phi <- log10(read.table(phi.file, sep=",", header=T)[, 2])
plot(NULL, NULL, ylim=range(expressionValues, na.rm = T) + c(-0.1, 0.1), xlim=range(obs.phi) + c(-0.1, 0.1),
     main = "Synthesis Rate", xlab = "true values", ylab = "estimated values")
upper.panel.plot(obs.phi[mixtureAssignment == 1], expressionValues[mixtureAssignment == 1], col="black")
legend("topleft", legend = paste("Mixture Element", 1:numMixtures),
       col = ribModel:::.mixtureColors[1:numMixtures], lty = rep(1, numMixtures), bty = "n")

#plot(parameter, what = "Mutation")
#plot(parameter, what = "Selection")
dev.off()


#plot(trace, what = "Expression", geneIndex = 999, mixture = 2)

#Plots the final aspects of the traces, that is, the mutation and selection traces
#for each amino acid and its set of synonymous codons.
mixture <- 1
pdf(paste(full.name, "CSP.pdf", sep=""), width = 11, height = 12)
plot(trace, what = "Mutation", mixture = 1)
plot(trace, what = "Selection", mixture = 1)
#plot(model, genome, parameter, samples = samples*0.1, mixture = 1, main = "Codon Usage Plot")

#TODO: Ask Jeremy about the rest of this code since it doesn't seem to do much in RStudio.
#However, the rest of the code before this is enough to judge how well the algorithm works.
names.aa <- aminoAcids()
selection.ci <- c()
mutation.ci <- c()
selection <- c()
mutation <- c()
codon.storage <- c()
csp.m <- read.table(mut.file, sep=",", header=T)
csp.e <- read.table(sel.file, sep=",", header=T)
#csp <- rbind(csp.m,csp.e)
#idx.eta <- 41:80
#idx.mu <- 1:40
ci <- .95
for(aa in names.aa)
{
  if(aa == "M" || aa == "W" || aa == "X") next
  codons <- AAToCodon(aa, T)
  codon.storage <- c(codon.storage,codons)
  for(i in 1:length(codons))
  {
    selection <- c(selection, parameter$getCodonSpecificPosteriorMean(mixture, samples*0.1, codons[i], 1, TRUE))
    sel.trace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 1)
    selection.ci <- c(selection.ci, quantile(tail(sel.trace, samples*0.1), probs = c(.5 - (ci / 2), .5 + (ci / 2)), names = FALSE))
    #selection.ci <- c(selection.ci, parameter$getCodonSpecificVariance(mixture, samples*0.1, codons[i], 1, TRUE))
    mutation <- c(mutation, parameter$getCodonSpecificPosteriorMean(mixture, samples*0.1, codons[i], 0, TRUE))
    mut.trace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 0)
    mutation.ci <- c(mutation.ci, quantile(tail(mut.trace, samples*0.1), probs = c(.5 - (ci / 2), .5 + (ci / 2)), names = FALSE))
    #mutation.ci <- c(mutation.ci, parameter$getCodonSpecificVariance(mixture, samples*0.1, codons[i], 0, TRUE))
  }
}
selection.ci <- matrix(selection.ci, ncol = 2, byrow = TRUE)
mutation.ci <- matrix(mutation.ci, ncol = 2, byrow= TRUE)
plot(NULL, NULL, xlim=range(csp.m[,3], na.rm = T), ylim=range(mutation),
     main = "Mutation", xlab = "True values", ylab = "Estimated values")
upper.panel.plot(x = csp.m[,3], y = mutation, sd.y = mutation.ci)
plot(NULL, NULL, xlim=range(csp.e[,3], na.rm = T), ylim=range(selection),
     main = "Selection", xlab = "True values", ylab = "Estimated values")
upper.panel.plot(x = csp.e[,3], y = selection, sd.y = selection.ci)
dev.off()
