#' config-rifampicin-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' Simulations to evaluate the MIC.fraction values for a range of kappa values.

#' The parameters defined here will be passed to functions in the
#' simulate-script.R file, which runs simulations either locally or on EULER
#' (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "rifampicin"
cell.type <- "WT"
output <- "survival"
total.sims <- 1000
sims.per.node <- 20
queue <- "24:00"
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    drug.outsideConcentration.MULT = 1,
    MIC.fraction = c(seq(0.28, 0.304, 0.002), seq(0.31, 0.35, 0.02)),
    kappa = c(0.5, 1.0, 2.0, 5.0, 10.0)
)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(drug, "survival", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "MIC-fraction/rifampicin-data/"
