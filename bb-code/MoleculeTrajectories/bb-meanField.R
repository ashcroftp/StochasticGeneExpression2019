#' bb-meanField.R
#' Author: Peter Ashcroft, ETH Zurich

#' A script to calculate some mean-field properties of cells

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/meanFieldModel.R")
source("R/dataAnalysisFunctions.R")
source("R/plotFunctions.R")

# Simulate ----
#' Load parameters & modify
drug <- "ciprofloxacin"
#drug <- "rifampicin"
#drug <- "cefotaxime"
basic.params <- loadParameters(drug)
basic.params$drug.outsideConcentration <- 1.0 * basic.params$MIC
#basic.params$rate.efflux.DNA.switchOn <- 10 * basic.params$rate.efflux.DNA.switchOn
#basic.params$rate.target.binding <- 0
#basic.params$rate.efflux.binding <- 0


#' Define some initial cell types and then create an explicit vector
cell.type <- "WT"
initial.cell <- c(1)
names(initial.cell) <- cell.type
#' Set the seed for reproducibility
set.seed(11)
#' Cell initial conditions
initial.cell.data <- simulateIC(basic.params, initial.cell)
#' Calculate ODE solution
t.max <- 400
ode.sol <- list(
    division = withinCellDF(generation = 0, data = integrateDivisionEqns(t.max, params = basic.params, IC = initial.cell.data), params = basic.params, add.volume = FALSE),
    dilution = withinCellDF(generation = 0, data = integrateDiluteEqns(t.max, params = basic.params, IC = initial.cell.data), params = basic.params, add.volume = FALSE)
)
#' Extract equilibrium value
ode.sol$dilution[nrow(ode.sol$dilution), "fraction.bound.targets"]
mean(ode.sol$division[ode.sol$division$time > t.max - 2 * basic.params$mean.generationTime, "fraction.bound.targets"])
# c(
#     min(ode.sol$division[ode.sol$division$time > t.max - 2 * basic.params$mean.generationTime, "fraction.bound.targets"]),
#     mean(ode.sol$division[ode.sol$division$time > t.max - 2 * basic.params$mean.generationTime, "fraction.bound.targets"]),
#     median(ode.sol$division[ode.sol$division$time > t.max - 2 * basic.params$mean.generationTime, "fraction.bound.targets"]),
#     max(ode.sol$division[ode.sol$division$time > t.max - 2 * basic.params$mean.generationTime, "fraction.bound.targets"])
# )

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
melted <- DFapply(names(ode.sol), function(type) {
    df <- melt(data = ode.sol[[type]], id.vars = c("time", "generation"), measure.vars = variables, variable.name = "molecule", value.name = "number")
    df$type <- factor(type, levels = names(ode.sol))
    return(df)
})
ggplot(melted, aes(x = time, y = number, colour = type)) +
    geom_step() +
    facet_wrap(~molecule, scales = "free_y", nrow = 4, labeller = label_panels) +
    theme_bw() +
    theme(legend.position = "None",
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "time (minutes)", y = NULL)
