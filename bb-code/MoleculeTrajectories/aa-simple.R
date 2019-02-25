#' aa-simple.R
#' Author: Peter Ashcroft, ETH Zurich

#' A script to plot simple molecule trajectories

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/meanFieldModel.R")
source("R/dataAnalysisFunctions.R")
source("R/plotFunctions.R")

# Simulate ----
#' Load parameters & modify
drug <- "ciprofloxacin"
basic.params <- loadParameters(drug)
basic.params$drug.outsideConcentration <- 0.5 * basic.params$MIC

#' Define some initial cell types and then create an explicit vector
cell.type <- "WT"
initial.cell <- c(1)
names(initial.cell) <- cell.type

#' Set the seed for reproducibility
set.seed(11)

#' Cell initial conditions
initial.cell.data <- simulateIC(basic.params, initial.cell)

#' Simulate (and keep track of a single cell)
t.max <- 200
out <- simulateCells(basic.params, initial.cell.data, cell.type, t.max, single.cell = TRUE)
#' Mean-field ODE solution 
ode.sol <- integrateDivisionEqns(t.max, params = out$cell.params[["1"]], IC = out$cell.data[["1"]][1,])

#' Add concentrations and fraction of bound targets to each cell's timeseries, and then combine into a single DF
ssa.sol <- DFapply(names(out$cell.data), function(id) {
    #' Isolate this cell's data and parameters
    this.cell.data <- out$cell.data[[id]]
    this.cell.params <- out$cell.params[[id]]
    #' Add volume, concentrations, and fraction bound targets
    withinCellDF(generation = id, data = this.cell.data, params = this.cell.params, add.volume = TRUE)
})
#' Arrange generation factors
ssa.sol$generation <- factor(ssa.sol$generation, levels = names(out$cell.data))
#' Add concentrations etc. (volume is already added)
ode.sol <- withinCellDF(generation = 0, data = ode.sol, params = out$cell.params[["1"]], add.volume = FALSE)

# Plot ----
variables <- c(
    `number of free drug molecules` = "drug",
    `number of unbound target proteins` = "target.protein",
    `number of bound target proteins` = "target.dimer",
    `fraction of bound target proteins` = "fraction.bound.targets",
    `efflux DNA state` = "efflux.DNA",
    `number of efflux mRNA` = "efflux.mRNA",
    `number of efflux proteins` = "efflux.protein",
    `number of bound efflux proteins` = "efflux.dimer"
)
plotWithinCell(ssa.data = ssa.sol, ode.data = ode.sol, measure.vars = variables) +
    facet_wrap(~molecule, scales = "free_y", nrow = 4, labeller = label_panels)
