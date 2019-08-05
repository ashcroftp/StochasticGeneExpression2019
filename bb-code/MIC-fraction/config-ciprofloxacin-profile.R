#' config-ciprofloxacin-profile.R
#' Author: Peter Ashcroft, ETH Zurich

#' Simulations to evaluate the survival probability profiles of WT cells for kappa-MIC.fraction pairs

#' The parameters defined here will be passed to functions in the
#' simulate-script.R file, which runs simulations either locally or on EULER
#' (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "survival"
total.sims <- 1000
sims.per.node <- 100
queue <- "04:00"
#' Variable parameters for which all combinations will be considered
var.params <- data.frame(
    kappa = c(1, 2, 3, 4, 5),
    MIC.fraction = c(0.08321406, 0.08058672, 0.07913400, 0.07814015, 0.07751760)
)
var.params <- lapply(seq(0,2,0.1), function(drug) cbind(var.params, data.frame(drug.outsideConcentration.MULT = drug)))
var.params <- do.call(rbind, var.params)
#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste(drug, "profile", sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "MIC-fraction/ciprofloxacin-data/"
