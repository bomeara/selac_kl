library(selac)
library(ape)

load("../data/yeast.practiceSELACunrestgamma106.Rdata")

start.rep = 1
stop.rep = 2
true.model = result
SimulateData(true.model=true.model, start.rep=start.rep, stop.rep=stop.rep)

for (sim.index in start.rep:stop.rep) {
  AnalyzeDataWithKnownValues(true.model, paste("fastaSet",sim.index, sep="_")
}
