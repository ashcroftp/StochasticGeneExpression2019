#' controlNoise.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("Simulate/filenames.R")

#' Type of drug and cells
drug <- "ciprofloxacin"
cell.type <- "STRUCT-BIND"

#' Directory where our files are stored
distribution.directory <- "MoleculeDistributions/STRUCT-data/"
survival.directory <- "PopulationSurvival/STRUCT-data/"
#' Load the protein distribution data
distribution.data <- read.csv(combined.file(distribution.directory, paste(drug, cell.type, sep = "_")))
distribution.data <- data.frame(
    reg.effect = factor(distribution.data$rate.efflux.DNA.switchOn.MULT),
    efflux.protein = distribution.data$efflux.protein
)
#' Load the survival data and extract the survival probability
survival.data <- read.csv(combined.file(survival.directory, paste(drug, cell.type, sep = "_")))
survival.data <- read.csv(paste0(survival.directory, paste(drug, "survival", sep = "_"), ".dat"))
survival.data <- data.frame(
    reg.effect = factor(survival.data$rate.efflux.DNA.switchOn.MULT),
    drug.dose = survival.data$drug.outsideConcentration,
    survive = survival.data$survive
)

#' Plots
protein.plot <- ggplot(distribution.data, aes(x = reg.effect, y = efflux.protein)) +
    geom_violin(aes(fill = reg.effect, colour = reg.effect), scale = "width", adjust = 2, position = position_dodge(0.9)) +
    scale_colour_manual(values = gradientPalette(length(levels(distribution.data$reg.effect)), palette = "Blues")) +
    scale_fill_manual(values = gradientPalette(length(levels(distribution.data$reg.effect)), palette = "Blues")) +
    stat_summary(fun.y = mean, geom = "point", fill = "black", shape = 95, size = 10, position = position_dodge(0.9)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None"
    ) +
    
    labs(x = "regulatory effect", y = "number of efflux protein per cell")

survival.plot <- ggplot(survival.data[survival.data$drug.dose <= 4, ], aes(x = drug.dose, y = survive, group = reg.effect)) +
    geom_line(aes(colour = reg.effect), size = 1.2) +
    geom_point(aes(colour = reg.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$reg.effect)), palette = "Blues"), name = "regulatory effect") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        legend.position = "None"
    ) +
    labs(x = "drug concentration (x MIC)", y = "survival probability")


# Combine plots ----
plot_grid(protein.plot, survival.plot, align = "b", axis = "b", nrow = 1, labels = c("AUTO"), label_size = 12, label_fontface = "plain", rel_widths = c(1,1.5))
