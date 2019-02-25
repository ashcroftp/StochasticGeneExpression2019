#' config-reg-STRUCT-BIND.R
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
cell.type <- "STRUCT-BIND"
output <- "survival"
total.sims <- 1000
sims.per.node <- 200
#' Variable parameters for which all combinations will be considered
effect <- c(1,2,8,32)
var.params <- expand.grid(
    drug.outsideConcentration.MULT = c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5)),
    mutant.effect = c(50),
    rate.efflux.DNA.switchOn.MULT = effect
)
var.params$rate.efflux.protein.translation.MULT <- 1/var.params$rate.efflux.DNA.switchOn.MULT
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(drug, cell.type, sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "PopulationSurvival/STRUCT-data/"
