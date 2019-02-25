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
sims.per.node <- 200
#' Variable parameters for which all combinations will be considered
var.params <- data.frame(
    kappa = c(1, 2, 3, 4, 5),
    MIC.fraction = c(0.10430000, 0.08951271, 0.08567227, 0.08319832, 0.08220000)
)
var.params <- lapply(c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5)), function(drug) cbind(var.params, data.frame(drug.outsideConcentration.MULT = drug)))
var.params <- do.call(rbind, var.params)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(drug, cell.type, sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "PopulationSurvival/kappa-data/"
