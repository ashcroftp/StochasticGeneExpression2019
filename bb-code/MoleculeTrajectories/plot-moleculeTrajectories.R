#' moleculeTrajectories.R
#' Author: Peter Ashcroft, ETH Zurich

#' A script to investigate within-cell dynamics in the absence and then presence of drug.
#' We can also compare with the mean-field ODE model.

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
basic.params$sd.generationTime <- 0
basic.params$rate.efflux.DNA.switchOn <- 2 * basic.params$rate.efflux.DNA.switchOn

#' Define some initial cell types and then create an explicit vector
cell.type <- "WT"
initial.cell <- c(1)
names(initial.cell) <- cell.type

#' Set the seed for reproducibility
seed <- 146
set.seed(seed)

#' Cell initial conditions
initial.cell.data <- simulateIC(basic.params, initial.cell)

#' Simulation without drug
basic.params$drug.outsideConcentration <- 0
#' Simulate (and keep track of a single cell)
t.max <- 136
out <- simulateCells(basic.params, initial.cell.data, cell.type, t.max, single.cell = TRUE)
#' Mean-field ODE solution 
ode.sol <- integrateDivisionEqns(t.max, params = out$cell.params[["1"]], IC = out$cell.data[["1"]][1,])

#' Simulation with drug
basic.params$drug.outsideConcentration <- basic.params$MIC
#' Identify last cell, which we will use to initialise the simulations
new.cell <- names(out$cell.data)[length(out$cell.data)]
new.cell.data <- out$cell.data[[new.cell]][1, ]
new.time <- new.cell.data$time
new.cell.data <- new.cell.data[, -1]
#' Simulate
new.out <- simulateCells(basic.params, new.cell.data, cell.type, t.max, single.cell = TRUE)
new.names <- as.character(as.numeric(names(new.out$cell.data)) + as.numeric(new.cell) - 1)
names(new.out$cell.data) <- new.names
names(new.out$cell.params) <- new.names
#' Mean-field ODE solution 
ode.ic <- ode.sol[ode.sol$time > new.time, ]
ode.ic <- ode.ic[1, ]
#ode.ic <- ode.sol[nrow(ode.ic), ]
new.ode.sol <- integrateDivisionEqns(t.max, params = new.out$cell.params[[new.cell]], IC = ode.ic)
#' Update time
for (id in new.names) new.out$cell.data[[id]]$time <- new.out$cell.data[[id]]$time + new.time
rm(id)
new.ode.sol$time <- new.ode.sol$time + new.time


#' Combine old (without drug) and new (with drug) outputs
#' Drop the last cell from the without drug DF
out$cell.data <- out$cell.data[names(out$cell.data) != new.cell]
out$cell.params <- out$cell.params[names(out$cell.params) != new.cell]
ode.sol <- ode.sol[ode.sol$time < new.time, ]
#' Add the new (with drug) data to the existing datafarmes
out$cell.data <- c(out$cell.data, new.out$cell.data)
out$cell.params <- c(out$cell.params, new.out$cell.params)
ode.sol <- rbind(ode.sol, new.ode.sol)

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

#' Save the data to file
data.directory <- "MoleculeTrajectories/zz-data/"
write.csv(ssa.sol, file = paste0(data.directory, "trajectory-SSA.dat"), row.names = F)
write.csv(ode.sol, file = paste0(data.directory, "trajectory-ODE.dat"), row.names = F)
list2csv(basic.params, paste0(data.directory, "trajectory.params"))
#' Save the working environment
sink(file =  paste0(data.directory, "trajectory.info"))
print(
    list(
        drug = drug,
        cell.type = cell.type,
        drug.on = new.time,
        drug.dose = basic.params$drug.outsideConcentration,
        timestamp = timestamp(quiet = T),
        version = R.Version()
    )
)
sink()

# Plot ----
#' Read data
data.directory <- "MoleculeTrajectories/zz-data/"
new.time <- 136
ssa.sol <- read.csv(paste0(data.directory, "trajectory-SSA.dat"))
ssa.sol$generation <- factor(ssa.sol$generation, levels = sort(unique(ssa.sol$generation)))
ode.sol <- read.csv(paste0(data.directory, "trajectory-ODE.dat"))
ode.sol$generation <- factor(ode.sol$generation, levels = sort(unique(ode.sol$generation)))

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
p <- plotWithinCell(ssa.data = ssa.sol, ode.data = ode.sol, measure.vars = variables) +
    # facet_wrap(~molecule, scales = "free_y", nrow = 4, dir = "v", labeller = label_panels)
    facet_wrap(~molecule, scales = "free_y", nrow = 4, labeller = label_panels)
#' Add shading for drug
p$layers <- c(
    geom_rect(data = data.frame(xmin = new.time, xmax = Inf, ymin = -Inf, ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.3, inherit.aes = F),
    p$layers
)
p
