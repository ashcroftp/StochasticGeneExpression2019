#' config-REG-ON.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "REG-ON"
output <- "molecules"
total.sims <- 10000
sims.per.node <- 1000
queue <- "04:00"

#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    drug.outsideConcentration = 0,
    mutant.effect = c(2,10,50,200)
)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(cell.type, "molecules", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "MoleculeDistributions/zz-data/"
