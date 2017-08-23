
#' @title Simulate data using known values
#'
#' @description
#' Creates a set of fasta files, one set per folder.
#'
#' @param true.model The result from a SelacOptimize() run
#' @param start.rep Index for first sim
#' @param end.rep Index for last sim
#' @param numgenes How many genes to simulate
SimulateData <- function(true.model, start.rep=1, end.rep=100, numgenes=5){

    phy = true.model$phy
    pars.mat <- true.model$mle.pars
    nuc.model = true.model$nuc.model

    if(numgenes == 5){
        rows.to.sample = 1:5
    }
    if(numgenes == 10){
        rows.to.sample = 1:10
    }
    if(numgenes == 20){
        rows.to.sample = order(pars.mat[,1])[seq(5,105,5)]
    }
    for(nrep.index in start.rep:end.rep){
        system(paste("mkdir", paste("fastaSet",nrep.index, sep="_"), sep=" "))
        for(gene.no.index in 1:numgenes){
            codon.freq.by.aa = true.model$codon.freq.by.aa[[rows.to.sample[gene.no.index]]]
            codon.freq.by.gene = true.model$codon.freq.by.gene[[rows.to.sample[gene.no.index]]]
            par.vec = c(pars.mat[rows.to.sample[gene.no.index],],1)
            optimal.aa = true.model$aa.optim[[rows.to.sample[gene.no.index]]]
            tmp <- SelacSimulator(phy=phy, pars=par.vec, aa.optim_array=optimal.aa, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, aa.properties=NULL, nuc.model=nuc.model, include.gamma=TRUE, k.levels=0, ncats=4, diploid=FALSE)
            write.dna(tmp, file=paste(paste("fastaSet",nrep.index, sep="_"), "/gene",  gene.no.index, ".fasta", sep=""), format="fasta", colw=1000000)
        }
    }
}

#' @title Analyze from simulated data
#'
#' @description
#' Loads and runs
#'
#' @param true.model The result from a SelacOptimize() run
#' @param data.dir Where the fasta files reside
#' @return A selac object
AnalyzeDataEstimatingValues <- function(true.model, data.dir) {
    result <- SelacOptimize(codon.data.path = data.dir, phy=true.model$phy, nuc.model=true.model$nuc.model)
    return(result)
}

#' @title Analyze from simulated data with known values
#'
#' @description
#' Loads and runs
#'
#' @param true.model The result from a SelacOptimize() run
#' @param data.dir Where the fasta files reside
#' @return A selac object
AnalyzeDataKnownValues <- function(true.model, data.dir) {
    result <- SelacOptimize(codon.data.path = data.dir, phy=true.model$phy, nuc.model=true.model$nuc.model)
    return(result)
}
