#' config-WT-growth.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "growth"
total.sims <- 100
sims.per.node <- 20
queue <- "24:00"
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    mutant.effect = c(1),
    drug.outsideConcentration.MULT = 10^seq(-1, 1, 0.1),
    efflux.binom.param = seq(0.5, 1, 0.1)
)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(cell.type, "growth", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "EffluxBias/zz-data/"
