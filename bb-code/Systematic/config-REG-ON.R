#' config-REG-ON.R
#' Author: Peter Ashcroft, ETH Zurich

#' The parameters defined here will be passed to functions in the simulate.R file,
#' which runs simulations either locally or on EULER (ETH HPC cluster).
#' Locations of config and data files are assigned at the bottom.

#' Required parameters
drug <- "ciprofloxacin"
cell.type <- "REG-ON"
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
        df <- params.output[c(2,3,5,6,8,9,12,14,15,17,18,24,27,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,101,102,104,105,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,128,129,131,132,134,135,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,276,279,285,288,297,325,326,327,328,329,330,331,332,335,336,337,338,339,340,341,345,347,348,353,354,355,356,357,358,359,360,362,363,364,365,366,367,368,369,372,374,375,377,378,406,407,408,409,410,411,412,413,415,416,417,418,419,420,421,422,424,425,426,427,428,429,430,431,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,549,558,567,598,599,600,601,602,608,609,610,611,626,627,629,630,635,636,638,639,644,645,647,648,649,650,651,652,653,654,655,658,659,660,661,662,663,664,668,669,670,671,672,676,677,678,679,680,681,682,683,685,686,687,688,689,690,691,692,695,696,697,698,699,700,701,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,722,723,724,725,726,727,728,729), ]
        df$drug.outsideConcentration.MULT <- drug
        return(df)
    }))
)

drug.values <- c(seq(60,100,10), seq(120,200,20))
var.params <- rbind(
    var.params,
    do.call(rbind, lapply(drug.values, function(drug) {
        df <- params.output[c(163,164,165,166,169,172,173,174,175,178,181,182,183), ]
        df$drug.outsideConcentration.MULT <- drug
        return(df)
    }))
)

drug.values <- seq(200,300,20)
var.params <- rbind(
    var.params,
    do.call(rbind, lapply(drug.values, function(drug) {
        df <- params.output[c(163,164,165,172,173,174,182,183), ]
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
