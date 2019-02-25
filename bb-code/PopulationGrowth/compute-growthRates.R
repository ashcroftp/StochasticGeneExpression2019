#' compute-growthRates.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("Simulate/filenames.R")

#' Drug and cell types to consider
drug <- "ciprofloxacin"
#drug <- "rifampicin"
#drug <- "cefotaxime"
cell.types <- c("REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
data.directory <- "PopulationGrowth/growth-data/"
#cell.types <- c("WT")
#data.directory <- "PopulationGrowth/kappa-data/"

#' Load the data and refactor
data <- DFapply(cell.types, function(cell.type) {
    df <- read.csv(combined.file(data.directory, paste(drug, cell.type, sep = "_")))
    df$cell.type <- factor(cell.type, levels = cell.types)
    return(df)
})
data$sim.label <- factor(data$sim.label, levels = sort(unique(data$sim.label)))
data$parameter.index <- factor(data$parameter.index, levels = sort(unique(data$parameter.index)))

#' Sample the data
sampled.data <- by(data, list(data$sim.label, data$parameter.index, data$cell.type), function(df) {
    DFapply(seq(0, max(data$sample.time), by = 10), function(t) {
        tmp <- df[df$sample.time <= t,]
        tmp <- tmp[tmp$sample.time == max(tmp$sample.time), ]
        tmp$sample.time <- t
        return(tmp)
    })
})
sampled.data <- do.call(rbind, sampled.data)

#' Sum over the simulations
summed.data <- by(sampled.data, list(sampled.data$sample.time, sampled.data$parameter.index, sampled.data$cell.type), function(df) {
    tmp <- df[df$sim.label == 1, ]
    tmp$num.cells <- sum(df$num.cells)
    return(tmp[, names(tmp)[names(tmp) != "sim.label"]])
})
summed.data <- do.call(rbind, summed.data)

#' Compute growth rates
growth.rates <- by(summed.data, list(summed.data$parameter.index, summed.data$cell.type), function(df) {
    df <- df[df$num.cells > 0, ]
    ic <- df[df$sample.time == 0, "num.cells"]
    fit.df <- data.frame(sample.time = df$sample.time, fit.var = log10(df$num.cells) - log10(ic))
    rate <- lm(fit.var ~ sample.time - 1, fit.df)$coefficients[["sample.time"]]

    tmp <- df[1, names(df)[!names(df) %in% c("sample.time", "num.cells")]]
    tmp$rate <- rate
    return(tmp)
})
growth.rates <- do.call(rbind, growth.rates)

#' Export data
write.csv(sampled.data, paste0(data.directory, paste(drug, "sampled", sep = "_"), ".dat"), row.names = F)
write.csv(summed.data, paste0(data.directory, paste(drug, "summed", sep = "_"), ".dat"), row.names = F)
write.csv(growth.rates, paste0(data.directory, paste(drug, "growth-rates", sep = "_"), ".dat"), row.names = F)
