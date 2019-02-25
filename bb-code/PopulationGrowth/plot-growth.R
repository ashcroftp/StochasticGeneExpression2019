#' plot-growth.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")

#' Drug and cell types to consider
drug <- "ciprofloxacin"
cell.types <- c("REG-ON", "REG-OFF", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
data.directory <- "PopulationGrowth/growth-data/"

#' Load the data
sampled.data <- read.csv(paste0(data.directory, paste(drug, "sampled", sep = "_"), ".dat"))
summed.data <- read.csv(paste0(data.directory, paste(drug, "summed", sep = "_"), ".dat"))
growth.rates <- read.csv(paste0(data.directory, paste(drug, "growth-rates", sep = "_"), ".dat"))



#' Plot time-kill curves
summed.data <- summed.data[summed.data$mutant.effect %in% c(1,200) & summed.data$drug.outsideConcentration.MULT %in% c(0.1, 1, 10, 100) & summed.data$cell.type %in% cell.types, ]
summed.data$cell.type <- factor(summed.data$cell.type, levels = cell.types)
summed.data$drug.outsideConcentration.MULT <- factor(summed.data$drug.outsideConcentration.MULT, levels = sort(unique(summed.data$drug.outsideConcentration.MULT)))

ggplot(summed.data, aes(x = sample.time, y = num.cells, colour = drug.outsideConcentration.MULT)) +
    geom_hline(data = summed.data[summed.data$sample.time == 0, ], aes(yintercept = num.cells), colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    #geom_line(aes(y = fit, colour = drug.dose)) +
    scale_color_manual(values = gradientPalette(length(levels(summed.data$drug.outsideConcentration.MULT)), palette = "Blues"), name = "drug concentration (x MIC)") +
    scale_y_log10(labels = sciFormat) +
    facet_wrap(~ cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    ) +
    guides(colour = guide_legend(nrow = 1)) +
    labs(x = "time (minutes)", y = "number of cells")



#' Plot growth rates
growth.rates <- growth.rates[growth.rates$cell.type %in% cell.types, ]
growth.rates$cell.type <- factor(growth.rates$cell.type, levels = cell.types)
growth.rates$mutant.effect <- factor(growth.rates$mutant.effect, levels = sort(unique(growth.rates$mutant.effect)))

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 0.75)
)
rrr.data <- with(
    c(rrr.params[[drug]], list(dose = 10^seq(log10(min(growth.rates$drug.outsideConcentration.MULT)), log10(max(growth.rates$drug.outsideConcentration.MULT)), length.out = 100))),
    data.frame(cell.type = "WT", mutant.effect = 1, drug.outsideConcentration.MULT = dose, rate =  psi.max - ((psi.max - psi.min) * dose^kappa)/(dose^kappa - psi.min/psi.max) )
)

ggplot(growth.rates, aes(x = drug.outsideConcentration.MULT, y = rate, colour = mutant.effect, group = mutant.effect)) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.min, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = rrr.data, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(growth.rates$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = sciFormat) +
    facet_wrap(~ cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "drug concentration (x MIC)", y = "log-10 growth rate (per minute)")
