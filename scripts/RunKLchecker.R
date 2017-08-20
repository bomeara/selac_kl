library(KLchecker)
library(phangorn)
library(ape)
library(geiger)
library(TreeSim)
rm(list=ls())


#
# simulator.true.model <- function(parameters) {
#   return(rnorm(n=1000, mean=parameters[1], sd=parameters[2]))
# }
#
# true.parameters <- c(1000, 0.12)
#
# param.estimator.candidate.model <- function(x) {
#   return(mean(1/x)) #here our candidate is exp
# }
#
# simulator.candidate.model <- function(parameters) {
#   return(rexp(n=1000, rate=parameters[1]))
# }
#
# likelihood.data.with.true.model <- function(data, parameters) {
#   return(sum(dnorm(data, mean=parameters[1], sd=parameters[2]), log=TRUE))
# }
#
# likelihood.data.with.candidate.model <- function(data) {
#   return(sum(dexp(data, rate=param.estimator.candidate.model(data)), log=TRUE))
# }
#
# K <- 1
#
# NormVsExp <- EstimateKL(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihoods.data.with.candidate.model, K)




simulator.true.model <- function(parameters) {
  sequences <- simSeq(parameters)
  final.object <- pml(tree=parameters$tree, data=sequences)
  return(final.object)
}


param.estimator.candidate.model <- function(x) {
  return(optim.pml(x, optNni=FALSE, model="JC", optGamma=FALSE, optQ=TRUE, optBf=TRUE, optRate = TRUE))
}

simulator.candidate.model <- function(fit) {
  return(simSeq(fit))
}

likelihood.data.with.true.model <- function(data, fit) {
  fit$data <- data
  return(optim.pml(fit, optNni = FALSE, optBf = FALSE, optQ = FALSE,
  optInv = FALSE, optGamma = FALSE, optEdge = FALSE, optRate = FALSE,
  optRooted = FALSE , max.it=0)$logLik)
}

likelihood.data.with.candidate.model <- function(data) {

  return(param.estimator.candidate.model(data)$logLik)
}

K <- 1

nchar.vector <- c(100, 1000, 3000, 5000)
ntax.vector <- c(25, 50, 100, 250)
results <- data.frame()
data(Laurasiatherian)
seed.tree = nj(dist.ml(Laurasiatherian))
seed.fit = pml(seed.tree, Laurasiatherian, k=1)
seed.fit = optim.pml(seed.fit, optNni=TRUE, model="HKY", optQ=TRUE)

for (rep in sequence(10)) {
  for (i in sequence(length(nchar.vector))) {
    for (j in sequence(length(ntax.vector))) {
      #
      #new.data <- Laurasiatherian
    #  new.data <- as.phyDat(as.character(new.data)[,1:nchar.vector[i]])
      #tree = nj(dist.ml(new.data))
    #  tree <- geiger::drop.random(tree, Ntip(tree) - ntax.vector[j])
      try({
        tree <- TreeSim::sim.bd.taxa(n=ntax.vector[j], numbsim=1, lambda=1, mu=0.9, frac=0.5, complete=FALSE)[[1]]
        tree$edge.lengths <- tree$edge.lengths/(max(branching.times(tree))) #so have less evolution
        new.data <- phangorn::simSeq(x=tree, l=nchar.vector[i], Q=seed.fit$Q, bf=seed.fit$bf, type="DNA" , rate=0.01)
        fit = pml(tree, new.data)
        fit = optim.pml(fit, optNni=TRUE, model="HKY", optGamma=FALSE, optQ=TRUE)
        true.parameters <- fit
        JCvsHKY <- EstimateKL(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihoods.data.with.candidate.model, K, nrep.outer=50, nrep.inner=50)
        results <- rbind(results, data.frame(nchar=nchar.vector[i], ntax=ntax.vector[j], KL=JCvsHKY[1], AIC=JCvsHKY[2], rep=rep))
        save(results, file="~/Dropbox/KL.rda")
        print(results)
      })
    }
  }
}
