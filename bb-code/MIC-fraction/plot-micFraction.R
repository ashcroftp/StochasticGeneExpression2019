#' plot-micFraction.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("Simulate/filenames.R")

#' Define drug and cell type
drug <- "ciprofloxacin"
#drug <- "rifampicin"
#drug <- "cefotaxime"
cell.type <- "WT"

#' Survival vs MIC fraction
#' Load data
data.directory <- "MIC-fraction/kappa-data/"
mic.data <- read.csv(paste0(data.directory, paste(drug, "survival", sep = "_"), ".dat"))
#' Refactor
#mic.data$kappa <- factor(format(mic.data$kappa, nsmall = 1))
mic.data$kappa <- factor(mic.data$kappa)
#' Plot
MIC.plot <- ggplot(mic.data, aes(x = MIC.fraction, y = survive, colour = kappa)) +
    geom_hline(yintercept = 0.5, colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(mic.data$kappa)), palette = "Blues"), name = parse(text = "kappa")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10)
    ) +
    labs(x = parse(text = "rho[MIC]"), y = "survival probability")

legend <- get_legend(MIC.plot)
MIC.plot <- MIC.plot + theme(legend.position = "none")

#' MIC fraction vs kappa
#' Interpolate the MIC-fraction value from this data
mic.fraction <- by(mic.data, list(mic.data$kappa), function(df) {
    sol <- approx(y = df$MIC.fraction, x = df$survive, xout = 0.5, ties = min)
    data.frame(kappa = unique(df$kappa), MIC.fraction = sol$y)
})
mic.fraction <- do.call(rbind, mic.fraction)
#mic.fraction$kappa <- as.numeric(as.character(mic.fraction$kappa))
#' Plot
kappa.plot <- ggplot(mic.fraction, aes(x = as.numeric(as.character(mic.fraction$kappa)), y = MIC.fraction, colour = kappa)) +
    geom_point(size = 2.0, shape = 16, show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(mic.fraction$kappa)), palette = "Blues"), name = parse(text = "kappa")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none"
    ) +
    labs(x = parse(text = "kappa"), y = parse(text = "rho[MIC]"))

plot_grid(
    plot_grid(
        MIC.plot + ggtitle("A: survival probability at 1x MIC"),
        kappa.plot + ggtitle("B: MIC fraction as a function of shape"),
        ncol = 2,
        align = "vh"
    ), legend, nrow = 1, rel_widths = c(2,0.2)
)




## Also plot growth rate and survival profiles ----

#' Growth rate vs kappa
data.directory <- "PopulationGrowth/kappa-data/"
#' Load the data
growth.rates <- read.csv(paste0(data.directory, paste(drug, "growth-rates", sep = "_"), ".dat"))
growth.rates$kappa <- factor(growth.rates$kappa)
growth.rates$drug.dose <- growth.rates$drug.outsideConcentration.MULT

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 2.5),
    cefotaxime = list(psi.max = 0.80/60, psi.min = -4.0/60, kappa = 0.75)
)
rrr.data <- with(
    c(rrr.params[[drug]], list(dose = 10^seq(log10(min(growth.rates$drug.outsideConcentration.MULT)), log10(max(growth.rates$drug.outsideConcentration.MULT)), length.out = 100))),
    data.frame(kappa = kappa, drug.dose = dose, rate =  psi.max - ((psi.max - psi.min) * dose^kappa)/(dose^kappa - psi.min/psi.max) )
)


growth.plot <- ggplot(growth.rates, aes(x = drug.dose, y = rate, colour = kappa, group = kappa)) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.min, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = rrr.data, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(growth.rates$kappa)), palette = "Blues"), name = parse(text = "kappa")) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = sciFormat) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none"
    ) +
    labs(x = "drug concentration (x MIC)", y = "base-10 growth rate (per minute)")


#' Survival curves
data.directory <- c("PopulationSurvival/kappa-data/")
#' Load the data
survival.data <- read.csv(paste0(data.directory, paste(drug, "survival", sep = "_"), ".dat"))
survival.data$kappa <- factor(survival.data$kappa)
#' Plot
survival.plot <- ggplot(survival.data[survival.data$drug.outsideConcentration.MULT <= 2, ], aes(x = drug.outsideConcentration.MULT, y = survive, colour = kappa)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$kappa)), palette = "Blues"), name = parse(text = "kappa")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none"
    ) +
    labs(x = "drug concentration (x MIC)", y = "survival probability")

#' Combined plots
plot_grid(
    plot_grid(
        MIC.plot + ggtitle("A: survival probability at 1x MIC"),
        kappa.plot + ggtitle("B: MIC fraction as a function of shape"),
        growth.plot + ggtitle("C: net population growth rate"),
        survival.plot + ggtitle("D: survival probability profile"),
        ncol = 2,
        align = "vh"
    ), legend, nrow = 1, rel_widths = c(2, 0.2)
)
