#' config-STRUCT-BIND-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "STRUCT-BIND"
output <- "survival"
total.sims <- 1000
sims.per.node <- 200
queue <- "04:00"
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    mutant.effect = c(2, 10, 50, 200),
    drug.outsideConcentration.MULT = c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5)),
    efflux.binom.param = seq(0.5, 1, 0.1)
)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(cell.type, "survival", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "EffluxBias/zz-data/"

