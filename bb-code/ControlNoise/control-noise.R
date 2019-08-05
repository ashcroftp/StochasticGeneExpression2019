#' control-noise.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Cell types and data directory
cell.type <- "STRUCT-BIND"
data.directory <- "ControlNoise/zz-data/"

## Load (and factor) the molecule distribution data ----
file.id <- "ControlNoise_molecules"
molecules.data <- read.csv(output.file(data.directory, file.id))
molecules.data <- data.frame(
    reg.effect = factor(molecules.data$rate.efflux.DNA.switchOn.MULT),
    efflux.protein = molecules.data$efflux.protein
)

#' Plot
molecules.plot <- ggplot(molecules.data, aes(x = reg.effect, y = efflux.protein)) +
    geom_violin(aes(fill = reg.effect, colour = reg.effect), scale = "width", adjust = 2, position = position_dodge(0.9)) +
    scale_colour_manual(values = gradientPalette(length(levels(molecules.data$reg.effect)), palette = "Greens")) +
    scale_fill_manual(values = gradientPalette(length(levels(molecules.data$reg.effect)), palette = "Greens")) +
    stat_summary(fun.y = mean, geom = "point", fill = "black", shape = 95, size = 10, position = position_dodge(0.9)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        plot.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None"
    ) +
    
    labs(x = "fold reduction of expression noise", y = "number of efflux protein per cell")



## Load (and factor) the lineage survival data ----
#' If the output file doesn't exist, we need to compute the survival probability
survival.file <- paste0(data.directory, "survival.dat")
if (!file.exists(survival.file) ) {
    file.id <- "ControlNoise_survival"
    survival.df <- computeSurvival(data.directory, file.id)
    #' Export
    write.csv(survival.df, survival.file, row.names = FALSE)
}

survival.data <- read.csv(survival.file)
survival.data <- data.frame(
    reg.effect = factor(survival.data$rate.efflux.DNA.switchOn.MULT),
    drug.dose = survival.data$drug.outsideConcentration.MULT,
    survive = survival.data$survive
)


survival.plot <- ggplot(survival.data[survival.data$drug.dose <= 5, ], aes(x = drug.dose, y = survive, group = reg.effect)) +
    geom_vline(xintercept = 1, colour = palette_OkabeIto[8]) +
    geom_line(aes(colour = reg.effect), size = 1.2) +
    geom_point(aes(colour = reg.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$reg.effect)), palette = "Greens"), name = "regulatory effect") +
    scale_x_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,10,1)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    coord_cartesian(xlim = c(0.4,5)) +
    theme_bw() +
    theme(
        plot.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "None"
    ) +
    labs(x = "drug concentration (x MIC)", y = "survival probability")

# Combine plots ----
plot_grid(
    molecules.plot + ggtitle("A: efflux protein distribution"),
    survival.plot + ggtitle("B: survival probability"),
    align = "b", axis = "b", nrow = 1, label_size = 12, label_fontface = "plain", rel_widths = c(1,1.1)
)
