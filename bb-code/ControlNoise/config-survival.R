#' config-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "STRUCT-BIND"
output <- "survival"
total.sims <- 1000
sims.per.node <- 100
queue <- "04:00"
#' Variable parameters for which all combinations will be considered
effect <- c(1,2,8,32)
var.params <- expand.grid(
    drug.outsideConcentration.MULT = c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5), seq(6,10,1)),
    mutant.effect = c(50),
    rate.efflux.DNA.switchOn.MULT = effect
)
var.params$rate.efflux.protein.translation.MULT <- 1/var.params$rate.efflux.DNA.switchOn.MULT
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- "ControlNoise_survival"
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "ControlNoise/zz-data/"
