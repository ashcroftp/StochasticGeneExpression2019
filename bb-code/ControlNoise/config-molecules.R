#' config-molecules.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "STRUCT-BIND"
output <- "molecules"
total.sims <- 10000
sims.per.node <- 1000
queue <- "04:00"

#' Variable parameters for which all combinations will be considered
effect <- c(1,2,8,32)
var.params <- expand.grid(
    drug.outsideConcentration.MULT = 0,
    mutant.effect = c(50),
    rate.efflux.DNA.switchOn.MULT = effect
)
var.params$rate.efflux.protein.translation.MULT <- 1/var.params$rate.efflux.DNA.switchOn.MULT
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- "ControlNoise_molecules"
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "ControlNoise/zz-data/"
