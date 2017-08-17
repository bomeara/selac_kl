library(KLchecker)
library(phangorn)
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
  return(optim.pml(x, optNni=TRUE, model="JC", optGamma=TRUE))
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

nchar.vector <- c(100, 1000, 3000)
results <- data.frame()
for (i in sequence(length(nchar.vector))) {
  data(Laurasiatherian)
  new.data <- Laurasiatherian
  new.data <- as.phyDat(as.character(new.data)[,1:nchar.vector[i]])
  tree = nj(dist.ml(new.data))
  fit = pml(tree, new.data, k=4)
  fit = optim.pml(fit, optNni=TRUE, model="HKY", optGamma=TRUE)
  true.parameters <- fit
  JCvsHKY <- EstimateKL(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihoods.data.with.candidate.model, K, nrep.outer=5, nrep.inner=4)
  results <- rbind(results, data.frame(nchar=nchar.vector[i], KL=JCvsHKY[1], AIC=JCvsHKY[2]))
  print(results)
}
