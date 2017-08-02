library(selac)
library(ape)

load("../data/yeast.practiceSELACunrestgamma106.Rdata")
selacSimControl(start.rep=76, end.rep=100, numgenes=10, yeast.fit.output=result)
