#' compute-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("Simulate/filenames.R")

#' Drug and cell types to consider
drug <- "ciprofloxacin"
#drug <- "rifampicin"
#drug <- "cefotaxime"
cell.types <- c("REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
data.directory <- "PopulationSurvival/survival-data/"
#cell.types <- c("WT")
#data.directory <- "PopulationSurvival/kappa-data/"
#cell.types <- c("STRUCT-BIND")
#data.directory <- "PopulationSurvival/STRUCT-data/"

#' Load the data
data <- DFapply(cell.types, function(cell.type) {
    df <- read.csv(combined.file(data.directory, paste0(drug, "_", cell.type)))
    df$cell.type <- factor(cell.type, levels = cell.types)
    return(df)
})
data$sim.label <- factor(data$sim.label, levels = sort(unique(data$sim.label)))
data$parameter.index <- factor(data$parameter.index, levels = sort(unique(data$parameter.index)))

#' Compute the survival probability and mean extinction time
survival <- by(data, list(data$parameter.index, data$cell.type), function(df) {
    tmp <- df[1, ]
    tmp$survive <- mean(df$survive)
    tmp$extinction.time <- mean(df[df$survive == 0, "exit.time"])
    tmp$extinction.time.sd <- sd(df[df$survive == 0, "exit.time"])
    return(tmp[, names(tmp)[!names(tmp) %in% c("sim.label", "num.cells", "exit.time")]])
})
survival <- do.call(rbind, survival)

#' Export
write.csv(survival, paste0(data.directory, paste(drug, "survival", sep = "_"), ".dat"), row.names = FALSE)
