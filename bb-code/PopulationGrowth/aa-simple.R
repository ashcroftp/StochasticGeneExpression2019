#' aa-simple.R
#' Author: Peter Ashcroft, ETH Zurich

#' Plot the growth rates of a population of cells over 1 hour

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/dataAnalysisFunctions.R")
source("R/plotFunctions.R")

#' Initialisation
#' Load parameters & modify
basic.params <- loadParameters(drug = "ciprofloxacin")
basic.params$drug.outsideConcentration <- basic.params$MIC

#' Define some initial cell types and then create an explicit vector
num.cells <- c(`WT` = 100)
initial.cell.types <- rep(names(num.cells), times = num.cells)

#' Cell initial conditions
initial.cell.data <- simulateIC(basic.params, num.cells)

#' Simulate
t.max <- 60
out <- simulateCells(basic.params, initial.cell.data, initial.cell.types, t.max)

#' Calculate population size and plot over time
pop.size <- analysePopSize(out$family.tree)
plotPopSize(pop.size)

