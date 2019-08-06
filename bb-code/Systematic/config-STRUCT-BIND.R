#' config-STRUCT-BIND.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "STRUCT-BIND"
output <- "survival"
total.sims <- 1000
sims.per.node <- 1000
queue <- "24:00"
#' We load the model parameters from Systematic/zz-data/params.output,
#' and then select the meaningful columns
params.output <- read.csv("./Systematic/zz-data/params.output")
params.output <- params.output[, c("rate.target.binding", "rate.efflux.binding", "rate.catalysis", "rate.drug.diffusion", "target.protein.total", "rate.target.protein.translation", "rate.efflux.mRNA.transcription", "MIC", "MIC.fraction")]
#' Add the mutant effect value(s) to the DF
params.output$mutant.effect <- c(200)
#' For each parameter set defined in the above DF, we consider multiple drug concentrations as multiples of the MIC
#' Note this column must be added to the DF after MIC
drug.values <- c(0, 0.5, 0.7, seq(0.9,1.4,0.1), seq(1.5,3,0.25), seq(3.5,5,0.5), seq(6,10,1))
var.params <- do.call(rbind, lapply(drug.values, function(drug) {
    df <- params.output
    df$drug.outsideConcentration.MULT <- drug
    return(df)
}))

drug.values <- c(seq(11,16,1), seq(18,30,2), seq(35,50,5))
var.params <- rbind(
    var.params,
    do.call(rbind, lapply(drug.values, function(drug) {
        df <- params.output[c(4,5,6,7,8,9,13,14,16,17,18,26,27,31,32,33,34,35,36,41,42,43,44,45,51,53,54,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,250,251,252,259,260,261,269,270,277,278,279,286,287,288,296,297,328,331,332,333,337,340,341,342,349,350,351,355,356,358,359,360,364,365,367,368,369,376,377,378,382,383,384,385,386,387,391,392,393,394,395,396,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,425,426,427,428,429,430,431,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,493,494,520,521,522,529,530,531,548,549,557,558,566,567,574,575,583,584,592,593,601,602,610,611,619,620,628,629,630,637,638,639,646,647,648,655,656,664,673,674,682,683,691,692,700,706,709,710,711,715,718,719,720,727,728,729), ]
        df$drug.outsideConcentration.MULT <- drug
        return(df)
    }))
)

drug.values <- c(seq(60,100,10), seq(120,200,20))
var.params <- rbind(
    var.params,
    do.call(rbind, lapply(drug.values, function(drug) {
        df <- params.output[c(8,9,17,18,27,36,45,85,88,89,90,94,97,98,99,107,108,112,113,114,115,116,117,121,122,123,124,125,126,131,132,134,135,144,152,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,251,278,279,287,288,297,331,332,340,341,350,358,359,360,367,368,369,377,378,386,387,395,396,405,412,413,421,422,430,431,436,439,440,441,448,449,450,457,458,459,463,464,465,466,467,468,472,473,474,475,476,477,481,482,483,484,485,486,549,558,574,601,602,610,628,629,630,637,638,639,647,648,655,664,682,691,709,710,718,719,727,728), ]
        df$drug.outsideConcentration.MULT <- drug
        return(df)
    }))
)

drug.values <- seq(220,300,20)
var.params <- rbind(
    var.params,
    do.call(rbind, lapply(drug.values, function(drug) {
        df <- params.output[c(89,90,98,99,116,117,125,126,135,169,170,178,179,188,196,197,198,205,206,207,215,216,223,224,225,230,233,234,240,243,359,412,448,468,476,477,485,486), ]
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
