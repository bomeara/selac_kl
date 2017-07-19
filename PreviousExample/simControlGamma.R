
library(selac)

selacSimControl <- function(start.rep=1, end.rep=100, phy, nuc.model, numcode=1, alpha=0.4, beta=0.1, Ne=5e6, numgenes=5, yeast.fit.output, aa.properties=NULL){
    
    phy = yeast.fit.output$phy
    pars.mat <- yeast.fit.output$mle.pars

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
            codon.freq.by.aa = yeast.fit.output$codon.freq.by.aa[[rows.to.sample[gene.no.index]]]
            codon.freq.by.gene = yeast.fit.output$codon.freq.by.gene[[rows.to.sample[gene.no.index]]]
            par.vec = c(pars.mat[rows.to.sample[gene.no.index],],1)
            optimal.aa = yeast.fit.output$aa.optim[[rows.to.sample[gene.no.index]]]
            tmp <- SelacSimulator(phy=phy, pars=par.vec, aa.optim_array=optimal.aa, codon.freq.by.aa=codon.freq.by.aa, codon.freq.by.gene=codon.freq.by.gene, numcode=1, aa.properties=NULL, nuc.model=nuc.model, include.gamma=TRUE, k.levels=0, ncats=4, diploid=FALSE)
            write.dna(tmp, file=paste(paste("fastaSet",nrep.index, sep="_"), "/gene",  gene.no.index, ".fasta", sep=""), format="fasta", colw=1000000)
        }
    }
}

load("yeast.practiceSELACunrestgamma106.Rdata")
#selacSimControl(start.rep=1, end.rep=25, nuc.model="UNREST", numgenes=10, yeast.fit.output=result)
#selacSimControl(start.rep=26, end.rep=50, nuc.model="UNREST", numgenes=10, yeast.fit.output=result)
#selacSimControl(start.rep=51, end.rep=75, nuc.model="UNREST", numgenes=10, yeast.fit.output=result)
selacSimControl(start.rep=76, end.rep=100, nuc.model="UNREST", numgenes=10, yeast.fit.output=result)



