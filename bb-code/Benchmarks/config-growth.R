#' config-growth.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Per config file, we only run one drug--cell.type combination (i.e. WT exposed to cipro).
#' Locations of config and data files are assigned at the bottom.


#' The following parameters must be specified in this file:
#'    - drug (one of c("ampicillin", "cefotaxime", "ciprofloxacin", "rifampicin", "tetracycline"))
#'    - cell.type (one of c("KO", "WT", "REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT"))
#'    - output (what type of output do we want from the simulations? One of c("survival))
#'    - total.sims (the number of simulations to perform per parameter set)
#'    - sims.per.node (how many simulations are sent to a given node? If zero simulations are run locally, otherwise this runs on EULER only)
#'    - var.params
#'    - file.id
#'    - data.directory
#' The parameters to be simulated over (var.params) are passed as a single dataframe.
#' The column names must be one of the parameter names, or be one of the names appended with .MULT.
#' Names with .MULT are multiplicative factors.
#' They are used, for example, to vary the drug concentration in units of MIC

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "growth"
total.sims <- 1
sims.per.node <- 0
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    drug.outsideConcentration.MULT = c(1)
)
#' File identifier (this is the base name, to which we append text and a suffix to dsitinguish the output files)
file.id <- paste(drug, cell.type, sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "Benchmarks/zz-data/"
