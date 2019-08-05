#' plot-micFraction.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Define drug and data directory
drug <- "ciprofloxacin"
#drug <- "rifampicin"
#drug <- "cefotaxime"
data.directory <- paste0("MIC-fraction/", drug, "-data/")

## Survival probability vs MIC fraction ----

#' If the output file doesn't exist, we need to compute the survival probability
survival.file <- paste0(data.directory, paste(drug, "survival.dat", sep = "_"))
if (!file.exists(survival.file) ) {
    file.id <- paste(drug, "survival", sep = "_")
    survival <- computeSurvival(data.directory, file.id)
    #' Export
    write.csv(survival, survival.file, row.names = FALSE)
}

#' Load data
mic.data <- read.csv(survival.file)
mic.data$kappa <- factor(mic.data$kappa)

#' Plot
MIC.plot <- ggplot(mic.data, aes(x = MIC.fraction, y = survive, colour = kappa)) +
    geom_hline(yintercept = c(0.1), colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(mic.data$kappa)), palette = "Reds"), name = parse(text = "kappa")) +
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
    sol <- fitSpline(x = df$MIC.fraction, y = df$survive, pred.y = 0.1)
    data.frame(kappa = unique(df$kappa), MIC.fraction = sol)
})
mic.fraction <- do.call(rbind, mic.fraction)

# spline.pred <- by(mic.data, list(mic.data$kappa), function(df) {
#     sol <- fitSpline(x = df$MIC.fraction, y = df$survive)
#     data.frame(kappa = unique(df$kappa), MIC.fraction = sol$x, survive = sol$y)
# })
# spline.pred <- do.call(rbind, spline.pred)
# MIC.plot + geom_line(data = spline.pred)

#' Plot
kappa.plot <- ggplot(mic.fraction, aes(x = as.numeric(as.character(mic.fraction$kappa)), y = MIC.fraction, colour = kappa)) +
    geom_point(size = 2.0, shape = 16, show.legend = F) +
    scale_x_continuous(limits = c(0,NA)) +
    scale_color_manual(values = gradientPalette(length(levels(mic.fraction$kappa)), palette = "Reds"), name = parse(text = "kappa")) +
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
        align = "vh",
        axis = "b"
    ), legend, nrow = 1, rel_widths = c(2,0.2)
)

## Growth rates vs kappa ----

#' If the output file doesn't exist, we need to compute the survival probability
growth.file <- paste0(data.directory, paste(drug, "growth.dat", sep = "_"))
if (!file.exists(growth.file) ) {
    file.id <- paste(drug, "growth", sep = "_")
    growth.rates <- computeGrowth(data.directory, file.id, final.times = c(60, 120, 180))
    #' Export
    write.csv(growth.rates, growth.file, row.names = FALSE)
}

#' Load the data
growth.rates <- read.csv(growth.file)
growth.rates$kappa <- factor(growth.rates$kappa)

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 2.5),
    cefotaxime = list(psi.max = 0.80/60, psi.min = -4.0/60, kappa = 0.75)
)
rrr.data <- with(
    c(rrr.params[[drug]], list(dose = 10^seq(log10(min(growth.rates$drug.outsideConcentration.MULT)), log10(max(growth.rates$drug.outsideConcentration.MULT)), length.out = 100))),
    data.frame(kappa = kappa, drug.outsideConcentration.MULT = dose, rate =  psi.max - ((psi.max - psi.min) * dose^kappa)/(dose^kappa - psi.min/psi.max) )
)


growth.plot <- ggplot(growth.rates[growth.rates$t.end == 180, ], aes(x = drug.outsideConcentration.MULT, y = rate, colour = kappa, group = kappa)) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.min, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    geom_line(data = rrr.data, linetype = "dashed", colour = palette_OkabeIto[8]) +
    scale_color_manual(values = gradientPalette(length(levels(growth.rates$kappa)), palette = "Reds"), name = parse(text = "kappa")) +
    scale_x_log10(breaks = c(0.1, 1, 10, 100), labels = sciFormat) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none"
    ) +
    labs(x = "drug concentration (x MIC)", y = "base-10 growth rate (per minute)")


## Survival profile for each kappa-MIC.fraction pair ----

#' If the output file doesn't exist, we need to compute the survival probability
profile.file <- paste0(data.directory, paste(drug, "profile.dat", sep = "_"))
if (!file.exists(profile.file) ) {
    file.id <- paste(drug, "profile", sep = "_")
    profile <- computeSurvival(data.directory, file.id)
    #' Export
    write.csv(profile, profile.file, row.names = FALSE)
}

#' Load data
profile.data <- read.csv(profile.file)
profile.data$kappa <- factor(profile.data$kappa)

#' Plot
profile.plot <- ggplot(profile.data[profile.data$drug.outsideConcentration.MULT <= 2, ], aes(x = drug.outsideConcentration.MULT, y = survive, colour = kappa)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(profile.data$kappa)), palette = "Reds"), name = parse(text = "kappa")) +
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
        profile.plot + ggtitle("D: survival probability profile"),
        ncol = 2,
        align = "vh"
    ), legend, nrow = 1, rel_widths = c(2, 0.2)
)


## Plot death rates as a function of kappa ----
death.rate <- DFapply(seq_len(nrow(mic.fraction)), function(i) {
    kappa <- as.numeric(levels(mic.fraction$kappa))[mic.fraction[i, "kappa"] ]
    MIC.fraction <- mic.fraction[i, "MIC.fraction"]
    rate.maxGrowth <- rrr.params[[drug]]$psi.max
    rate.minGrowth <- rrr.params[[drug]]$psi.min
    
    death.pars <- c(
        A = rate.maxGrowth * (rate.maxGrowth - rate.minGrowth) * (((1/MIC.fraction)^kappa) - 1) / (rate.maxGrowth * ((1/MIC.fraction)^kappa) - (rate.maxGrowth - rate.minGrowth)),
        B = rate.minGrowth * ((1/MIC.fraction)^kappa) / (rate.maxGrowth * ((1/MIC.fraction)^kappa) - (rate.maxGrowth - rate.minGrowth))
    )
    
    bound.fraction <- seq(0, 1, length.out = 101)
    death.rate <- (death.pars[["A"]] * ((bound.fraction / MIC.fraction)^kappa)) / (((bound.fraction / MIC.fraction)^kappa) - death.pars[["B"]])
    
    data.frame(
        kappa = mic.fraction[i, "kappa"],
        MIC.fraction = MIC.fraction,
        rate.maxGrowth = rate.maxGrowth,
        rate.minGrowth = rate.minGrowth,
        bound.fraction = bound.fraction,
        death.rate = death.rate,
        growth.rate = rate.maxGrowth - death.rate
    )
})

ggplot(death.rate, aes(x = bound.fraction, y = growth.rate, colour = kappa)) +
    geom_hline(aes(yintercept = rate.maxGrowth), colour = palette_OkabeIto[8]) +
    geom_hline(aes(yintercept = rate.minGrowth), colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) + 
    scale_color_manual(values = gradientPalette(length(levels(death.rate$kappa)), palette = "Reds"), name = parse(text = "kappa")) +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "top"
    ) +
    labs(x = "fraction of bound targets", y = "growth rate (per minute)")
