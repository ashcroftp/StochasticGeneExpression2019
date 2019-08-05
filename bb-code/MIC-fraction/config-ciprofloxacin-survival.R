#' config-ciprofloxacin-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' Simulations to evaluate the MIC.fraction values for a range of kappa values.

#' The parameters defined here will be passed to functions in the
#' simulate-script.R file, which runs simulations either locally or on EULER
#' (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "survival"
total.sims <- 1000
sims.per.node <- 100
queue <- "04:00"
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    drug.outsideConcentration.MULT = 1,
    MIC.fraction = c(seq(0.07, 0.085, 0.0005), seq(0.090, 0.120, 0.005)),
    kappa = seq(1,5)
)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(drug, "survival", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "MIC-fraction/ciprofloxacin-data/"
