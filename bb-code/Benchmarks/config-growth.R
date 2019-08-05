#' config-growth.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "KO"
output <- "growth"
total.sims <- 1
sims.per.node <- 0
queue <- "24:00"
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    drug.outsideConcentration.MULT = c(1)
)
#' File identifier (this is the base name, to which we append text and a suffix to dsitinguish the output files)
file.id <- paste(drug, cell.type, "growth", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "Benchmarks/zz-data/"
