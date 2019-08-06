#' config-WT.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "WT"
output <- "survival"
total.sims <- 1000
sims.per.node <- 1000
queue <- "24:00"
#' We load the model parameters from Systematic/zz-data/params.output,
#' and then select the meaningful columns
params.output <- read.csv("./Systematic/zz-data/params.output")
params.output <- params.output[, c("rate.target.binding", "rate.efflux.binding", "rate.catalysis", "rate.drug.diffusion", "target.protein.total", "rate.target.protein.translation", "rate.efflux.mRNA.transcription", "MIC", "MIC.fraction")]
#' Add the mutant effect value(s) to the DF
params.output$mutant.effect <- c(1)
#' For each parameter set defined in the above DF, we consider multiple drug concentrations as multiples of the MIC
#' Note this column must be added to the DF after MIC
drug.values <- c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5))
var.params <- do.call(rbind, lapply(drug.values, function(drug) {
    df <- params.output
    df$drug.outsideConcentration.MULT <- drug
    return(df)
}))

drug.values <- seq(5.5,15.0,0.5)
var.params <- rbind(
    var.params,
    do.call(rbind, lapply(drug.values, function(drug) {
        df <- params.output[c(163,164,165,172,173,174,181,182,183), ]
        df$drug.outsideConcentration.MULT <- drug
        return(df)
    }))
)

#' Set kappa value as constant
var.params$kappa <- 3.0

#' File identifier (this is the base name, to which we append text and a suffix to distinguish the output files)
file.id <- paste("systematic", cell.type, sep = "_")
#' Location of the data directory (relative to the "home directory" bb-code/)
data.directory <- "Systematic/zz-data/"
