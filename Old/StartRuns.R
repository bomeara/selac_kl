library(selac)
library(ape)
source("BasicSimulationScripts.R")

load("../data/yeast.practiceSELACunrestgamma106.Rdata")

start.rep = 1
end.rep = 20
true.model = result
rm(result)
SimulateData(true.model=true.model, start.rep=start.rep, end.rep=end.rep)

for (sim.index in start.rep:end.rep) {
  print(Sys.time())
  print(paste0("Starting run ", sim.index))
  result <- AnalyzeDataEstimatingValues(true.model, paste0(paste("fastaSet",sim.index, sep="_"),'/'))
  save(result, file=paste0("SimResult_",sim.index,".rda"))
}
