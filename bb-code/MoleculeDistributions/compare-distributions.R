#' compare-distributions.R
#' Author: Peter Ashcroft, ETH Zurich

#' A script to investigate distributions of protein molecules across cell types

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("Simulate/filenames.R")


#' Type of drug and cells
drug <- "ciprofloxacin"
cell.types <- c("WT", "REG-ON", "REG-OFF")

#' Load the distribution data
data.directory <- c(`WT` = "MoleculeDistributions/WT-data/", `REG-ON` = "MoleculeDistributions/REG-data/", `REG-OFF` = "MoleculeDistributions/REG-data/")
data <- DFapply(cell.types, function(cell.type) {
    data <- read.csv(combined.file(data.directory[[cell.type]], paste(drug, cell.type, sep = "_")))
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

cols <- palette_OkabeIto[c(7,5,6)]
names(cols) <- levels(data$cell.type)

dodge <- position_dodge(1.0)
ggplot(data, aes(x = mutant.effect, y = efflux.protein, fill = cell.type)) +
    geom_violin(aes(colour = cell.type), scale = "width", position = dodge) +
    scale_colour_manual(values = cols, name = "cell type") +
    scale_fill_manual(values = cols, name = "cell type") +
    stat_summary(aes(group = cell.type), fun.y = mean, geom = "point", fill = "black", shape = 95, size = 10, position = dodge) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "mutation effect", y = "number of efflux protein per cell")
