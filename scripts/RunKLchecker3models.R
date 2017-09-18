library(KLchecker)
library(phangorn)
library(ape)
library(geiger)
library(TreeSim)
rm(list=ls())
#options(error = utils::recover)


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
# NormVsExp <- EstimateKL(simulator.true.model, true.parameters, param.estimator.candidate.model, simulator.candidate.model, likelihood.data.with.true.model, likelihood.data.with.candidate.model, K)




simulator.true.model <- function(parameters) {
  sequences <- simSeq(parameters)
  final.object <- pml(tree=parameters$tree, data=sequences)
  return(final.object)
}


param.estimator.candidate.model.JC <- function(x) {
  return(optim.pml(x, optNni=FALSE, model="JC", optGamma=FALSE, optQ=TRUE, optBf=TRUE, optRate = TRUE, optRooted=TRUE))
}

param.estimator.candidate.model.HKY <- function(x) {
  return(optim.pml(x, optNni=FALSE, model="HKY", optGamma=FALSE, optQ=TRUE, optBf=TRUE, optRate = TRUE, optRooted=TRUE))
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

likelihood.data.with.candidate.model.HKY <- function(data) {
  return(param.estimator.candidate.model.HKY(data)$logLik)
}


likelihood.data.with.candidate.model.JC <- function(data) {
  return(param.estimator.candidate.model.JC(data)$logLik)
}


RunIndividual <- function(model.K.vector, k.index, param.estimator.candidate.model.JC, param.estimator.candidate.model.HKY, TreeSim, sim.bd.taxa, ntax.vector, j, phangorn, simSeq, nchar.vector, i, seed.fit, simulator.true.model, simulator.candidate.model, likelihood.data.with.true.model, likelihood.data.with.candidate.model.JC,  likelihood.data.with.candidate.model.HKY, results, rep) {
    K <- model.K.vector[k.index]
    param.estimator.candidate.model <- param.estimator.candidate.model.JC
    likelihood.data.with.candidate.model <- likelihood.data.with.candidate.model.JC
    if(k.index==2) {
        param.estimator.candidate.model <- param.estimator.candidate.model.HKY
        likelihood.data.with.candidate.model <- likelihood.data.with.candidate.model.HKY
    }
    new.data <- NULL
    tree <- NULL
    nrow.new.data <- 0
    attempts <- 1
    #    while(nrow.new.data<=3 & attempts<50) { #need variation
    tree <- TreeSim::sim.bd.taxa(n=ntax.vector[j], numbsim=1, lambda=1, mu=0.9, frac=0.5, complete=FALSE)[[1]]
    tree$edge.length <-  tree$edge.length/(max(branching.times(tree))) #so have less evolution; if too little, starts adding more
    new.data <- phangorn::simSeq(x=tree, l=nchar.vector[i], Q=seed.fit$Q, bf=seed.fit$bf, type="DNA" , rate=0.01 * attempts)
    nrow.new.data <- nrow(as.character(unique(new.data)))
    attempts <- attempts+1
    #print(attempts)
    #      if(attempts>=50) {
    #        new.data <- NULL
    #      }
    #    }
    #    if(!is.null(new.data)) {
    fit = pml(tree, new.data)
    
    fit = optim.pml(fit, optNni=FALSE, model="GTR", optGamma=FALSE, optQ=TRUE)
    true.parameters <- fit
    model.comparison <- EstimateKL(simulator.true.model=simulator.true.model, true.parameters=true.parameters, param.estimator.candidate.model=param.estimator.candidate.model, simulator.candidate.model=simulator.candidate.model, likelihood.data.with.true.model=likelihood.data.with.true.model, likelihoods.data.with.candidate.model=likelihood.data.with.candidate.model, K=K, nrep.outer=50, nrep.inner=50)
    results <- rbind(results, data.frame(nchar=nchar.vector[i], ntax=ntax.vector[j], KL=model.comparison[1], AIC=model.comparison[2], rep=rep, K=K))
    save(results, file="KLrun3model.rda")
    print(results)
    return(results)
}

nchar.vector <- c(3,4,5,10,20,50,100,200)
ntax.vector <- c(4,5,6,7,8,9,10, 25, 50, 100)
results <- data.frame()
data(Laurasiatherian)
model.K.vector <- c(1,2)
seed.tree = nj(dist.ml(Laurasiatherian))
seed.fit = pml(seed.tree, Laurasiatherian, k=1)
seed.fit = optim.pml(seed.fit, optNni=TRUE, model="GTR", optQ=TRUE)

for (rep in sequence(10)) {
  for (i in sequence(length(nchar.vector))) {
    for (j in sequence(length(ntax.vector))) {
      for (k.index in sequence(model.K.vector)) {

        #
        #new.data <- Laurasiatherian
      #  new.data <- as.phyDat(as.character(new.data)[,1:nchar.vector[i]])
        #tree = nj(dist.ml(new.data))
      #  tree <- geiger::drop.random(tree, Ntip(tree) - ntax.vector[j])
          results.this <- NULL
        try(results.this <- RunIndividual(model.K.vector=model.K.vector, k.index=k.index, param.estimator.candidate.model.JC=param.estimator.candidate.model.JC, param.estimator.candidate.model.HKY=param.estimator.candidate.model.HKY, TreeSim=TreeSim, sim.bd.taxa=sim.bd.taxa, ntax.vector=ntax.vector, j=j, phangorn=phangorn, simSeq=simSeq, nchar.vector=nchar.vector, i=i, seed.fit=seed.fit, simulator.true.model=simulator.true.model, simulator.candidate.model=simulator.candidate.model, likelihood.data.with.true.model=likelihood.data.with.true.model, likelihood.data.with.candidate.model.JC=likelihood.data.with.candidate.model.JC, likelihood.data.with.candidate.model.HKY=likelihood.data.with.candidate.model.HKY, results=results, rep=rep))
          if(!is.null(results.this)) {
            results <- results.this   
          }

      #    } else {
      #      print(paste0("Tried but failed nchar = ", nchar.vector[i], " ntax = ", ntax.vector[j], " model K = ", K))
    #      }
      }
    }
  }
}
