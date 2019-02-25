#' plot-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")

#' Drug and cell types to consider
drug <- "ciprofloxacin"
cell.types <- c("REG-ON", "REG-OFF", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
data.directory <- "PopulationSurvival/survival-data/"

#' Load the data
data <- read.csv(paste0(data.directory, paste(drug, "survival", sep = "_"), ".dat"))
#' Refactor
data <- data[data$cell.type %in% cell.types, ]
data$cell.type <- factor(data$cell.type, levels = cell.types)
data$parameter.index <- factor(data$parameter.index, levels = sort(unique(data$parameter.index)))
data$mutant.effect <- factor(data$mutant.effect, levels = sort(unique(data$mutant.effect)))
data$drug.off <- factor(data$drug.off, levels = sort(unique(data$drug.off)))


## Survival probabilities and extinction times ----
survival.data <- data[data$drug.off == Inf, ]

#' Plot the survival probabilities
ggplot(survival.data, aes(x = drug.outsideConcentration.MULT, y = survive, group = mutant.effect)) +
    geom_line(aes(colour = mutant.effect), size = 1.2) +
    geom_point(aes(colour = mutant.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "drug concentration (x MIC)", y = "survival probability")

IC <- by(survival.data, list(survival.data$mutant.effect, survival.data$cell.type), function(df) {
    sol <- approx(y = df$drug.outsideConcentration.MULT, x = df$survive, xout = c(0.5, 0.01, 0.001), ties = min)
    if (is.null(sol)) return(NULL)
    data.frame(
        cell.type = unique(df$cell.type),
        mutant.effect = unique(df$mutant.effect),
        IC50 = sol$y[[1]],
        IC99 = sol$y[[2]],
        IC99.9 = sol$y[[3]]
    )
})
do.call(rbind, IC)

#' Plot the extinction times (with standard deviation)
ggplot(survival.data, aes(x = drug.outsideConcentration.MULT, y = extinction.time, group = mutant.effect)) +
    geom_line(aes(colour = mutant.effect), size = 1.2) +
    geom_errorbar(aes(colour = mutant.effect, ymin = extinction.time - extinction.time.sd, ymax = extinction.time + extinction.time.sd), size = 0.1) +
    geom_point(aes(colour = mutant.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "drug concentration (x MIC)", y = "mean lineage extinction time")


## Pusled treatments ----
pulsed.data <- data[data$drug.outsideConcentration.MULT %in% c(1:5) & data$mutant.effect %in% c(1, 200), ]
pulsed.data$drug.outsideConcentration.MULT <- factor(pulsed.data$drug.outsideConcentration.MULT)
levels(pulsed.data$drug.off)[levels(pulsed.data$drug.off) == "Inf"] <- "1200"
#' Plot the survival probability across drug pulses
ggplot(pulsed.data, aes(x = drug.outsideConcentration.MULT, y = drug.off)) +
    geom_tile(aes(fill = survive)) +
    scale_fill_distiller(palette = "Blues", name = "survival probability", type = 'seq', na.value = "grey60", direction = 1, limits = c(0,1)) +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    labs(x = "drug concentration (x MIC)", y = "treatment duration (minutes)") +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0))
