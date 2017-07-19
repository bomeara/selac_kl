library(selac)
library(ape)

phy <- ape::stree(8, type="balanced")
phy$edge.length <- phy$edge.length/50 

simulations <- list()
nreps = 50


for (i in sequence(nreps)) {
  sim.result <- selac::SelacSimulator(phy, STUFF)
}

