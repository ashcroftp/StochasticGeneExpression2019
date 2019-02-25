#' config-cipro-kappa.R
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
output <- "survival"
total.sims <- 1000
sims.per.node <- 100
#' Variable parameters for which all combinations will be considered
var.params <- expand.grid(
    drug.outsideConcentration.MULT = 1,
    MIC.fraction = seq(0.07, 0.12, by = 0.0025),
    kappa = c(1:5)
)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(drug, cell.type, sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "MIC-fraction/kappa-data/"
