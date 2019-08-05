#' compare-WT-REG.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("Simulate/filenames.R")


#' Cell types and data directory
cell.types <- c("WT", "REG-ON", "REG-OFF")
data.directory <- "MoleculeDistributions/zz-data/"

#' Load the distribution data
data <- DFapply(cell.types, function(cell.type) {
    file.id <- paste(cell.type, "molecules", sep = "_")
    data <- read.csv(output.file(data.directory, file.id))
    df <- data.frame(
        cell.type = factor(cell.type, levels = cell.types),
        mutant.effect = factor(data$mutant.effect),
        efflux.protein = data$efflux.protein
    )
    return(df)
})

#' Refactor for wildtype label
levels(data$mutant.effect)[levels(data$mutant.effect) == "1"] <- "WT"
levels(data$cell.type)[levels(data$cell.type) == "WT"] <- "WT or STRUCT*"

cols <- cellColours()
names(cols)[names(cols) == "WT"] <- "WT or STRUCT*"

dodge <- position_dodge(1.0)
ggplot(data, aes(x = mutant.effect, y = efflux.protein, fill = cell.type)) +
    geom_violin(aes(colour = cell.type), scale = "width", position = dodge) +
    scale_colour_manual(values = cols, name = "cell type") +
    scale_fill_manual(values = cols, name = "cell type") +
    stat_summary(aes(group = cell.type), fun.y = mean, geom = "point", fill = "black", shape = 95, size = 10, position = dodge) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "mutation effect", y = "number of efflux protein per cell")
