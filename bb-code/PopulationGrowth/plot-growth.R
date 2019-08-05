#' plot-growth.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Data directory and cell types to consider
data.directory <- "PopulationGrowth/zz-data/"
cell.types <- c("REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
drug <- "ciprofloxacin"

## Read summed data (or generate the file first if missing), and plot time-kill curves ----
summed.file <- paste0(data.directory, paste(drug, "summed.dat", sep = "_"))
if (!file.exists(summed.file) ) {
    summed <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(drug, cell.type, "growth", sep = "_")
        sampled.df <- read.csv(output.file(data.directory, file.id))
        #' Add cell type to DF, and refactor
        sampled.df$cell.type <- factor(cell.type, levels = cell.types)
        sampled.df$sim.label <- factor(sampled.df$sim.label, levels = sort(unique(sampled.df$sim.label)))
        sampled.df$param.index <- factor(sampled.df$param.index, levels = sort(unique(sampled.df$param.index)))
        
        #' Now sum over the replicates (sim.label)
        summed.df <- by(sampled.df, list(sampled.df$sample.time, sampled.df$param.index), function(df) {
            tmp <- df[1, ]
            tmp$num.cells <- sum(df$num.cells)
            return(tmp[, names(tmp)[names(tmp) != "sim.label"]])
        })
        summed.df <- do.call(rbind, summed.df)
        return(summed.df)
    })
    #' Export
    write.csv(summed, summed.file, row.names = FALSE)
}
summed.data <- read.csv(summed.file)
summed.data$cell.type <- factor(summed.data$cell.type, levels = cell.types)
summed.data$drug.outsideConcentration.MULT <- factor(summed.data$drug.outsideConcentration.MULT, levels = sort(unique(summed.data$drug.outsideConcentration.MULT)))
#summed.data$drug <- log10(as.numeric(levels(summed.data$drug.outsideConcentration.MULT))[summed.data$drug.outsideConcentration.MULT])

#' Plot the time-kill curves for 1 mutant level and a subset of drug concentrations
summed.data <- summed.data[summed.data$mutant.effect %in% c(1,200) & summed.data$drug.outsideConcentration.MULT %in% c(0.1, 1, 10, 100) & summed.data$cell.type %in% cell.types, ]
summed.data$mutant.effect <- factor(summed.data$mutant.effect, levels = c(1,200))
summed.data$drug.outsideConcentration.MULT <- factor(summed.data$drug.outsideConcentration.MULT, levels = c(0.1, 1, 10, 100))

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 0.75)
)
max.growth <- with(
    c(rrr.params[[drug]], list(time = seq(0, max(summed.data$sample.time), length.out = 101))),
    DFapply(cell.types, function(cell.type) data.frame(cell.type = cell.type, mutant.effect = 1, drug.outsideConcentration.MULT = min(levels(summed.data$drug.outsideConcentration.MULT)), drug = 0, sample.time = time, num.cells =  10^5 * 10 ^ (psi.max * time), t.end = factor(60) ) )
)
max.growth.90 <- with(
    c(rrr.params[[drug]], list(time = seq(0, max(summed.data$sample.time), length.out = 101))),
    DFapply(cell.types, function(cell.type) data.frame(cell.type = cell.type, mutant.effect = 1, drug.outsideConcentration.MULT = min(levels(summed.data$drug.outsideConcentration.MULT)), drug = 0, sample.time = time, num.cells =  10^5 * 10 ^ (0.1 * psi.max * time), t.end = factor(60) ) )
)

ggplot(summed.data, aes(x = sample.time, y = num.cells, colour = drug.outsideConcentration.MULT)) +
    geom_hline(data = summed.data[summed.data$sample.time == 0, ], aes(yintercept = num.cells), colour = palette_OkabeIto[8]) +
    #geom_vline(data = summed.data[summed.data$sample.time %in% seq(0, max(summed.data$sample.time), 60), ], aes(xintercept = sample.time), colour = palette_OkabeIto[8]) +
    geom_line(data = max.growth, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = max.growth.90, linetype = "dashed", colour = palette_OkabeIto[8]) +
    
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    #scale_colour_viridis_c(limits = c(0, 1), name = "drug concentration (log10)") +
    #geom_line(aes(y = fit, colour = drug.dose)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(11, "PiYG")[5:1])(length(levels(summed.data$drug.outsideConcentration.MULT))), name = "drug concentration (x MIC)") +
    #scale_color_manual(values = gradientPalette(length(levels(summed.data$drug.outsideConcentration.MULT)), palette = "Greys"), name = "drug concentration (x MIC)") +
    scale_y_log10(labels = sciFormat) +
    facet_wrap(~ cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    #guides(colour = guide_legend(nrow = 1)) +
    labs(x = "time (minutes)", y = "number of cells")


# ggplot(summed.data, aes(x = sample.time, y = num.cells, colour = drug.outsideConcentration.MULT)) +
#     geom_hline(data = summed.data[summed.data$sample.time == 0, ], aes(yintercept = num.cells), colour = palette_OkabeIto[8]) +
#     geom_line(size = 1.2) +
#     geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
#     #geom_line(aes(y = fit, colour = drug.dose)) +
#     scale_color_manual(values = gradientPalette(length(levels(summed.data$drug.outsideConcentration.MULT)), palette = "Blues"), name = "drug concentration (x MIC)") +
#     scale_y_log10(labels = sciFormat) +
#     facet_wrap(~ cell.type, ncol = 2, labeller = label_panels) +
#     theme_bw() +
#     theme(
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "top"
#     ) +
#     guides(colour = guide_legend(nrow = 1)) +
#     labs(x = "time (minutes)", y = "number of cells")


## Compute growth rates from the summed data and plot dose-response curves ----
growth.file <- paste0(data.directory, paste(drug, "growth.dat", sep = "_"))
if (!file.exists(growth.file) ) {
    growth.df <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(drug, cell.type, "growth", sep = "_")
        growth.rates <- computeGrowth(data.directory, file.id, final.times = c(60, 120, 180))
        growth.rates$cell.type <- factor(cell.type, levels = cell.types)
        return(growth.rates)
    })
    #' Export
    write.csv(growth.df, growth.file, row.names = FALSE)
}
growth.data <- read.csv(growth.file)
growth.data <- growth.data[growth.data$cell.type %in% cell.types & growth.data$t.end == 180, ]
growth.data$cell.type <- factor(growth.data$cell.type, levels = cell.types)
growth.data$mutant.effect <- factor(growth.data$mutant.effect, levels = sort(unique(growth.data$mutant.effect)))
#growth.data$t.end <- factor(growth.data$t.end, levels = sort(unique(growth.data$t.end)))

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 0.75)
)
rrr.data <- with(
    c(rrr.params[[drug]], list(dose = 10^seq(log10(min(growth.data$drug.outsideConcentration.MULT)), log10(max(growth.data$drug.outsideConcentration.MULT)), length.out = 100))),
    data.frame(cell.type = "WT", mutant.effect = 1, drug.outsideConcentration.MULT = dose, rate =  psi.max - ((psi.max - psi.min) * dose^kappa)/(dose^kappa - psi.min/psi.max) )
)

ggplot(growth.data, aes(x = drug.outsideConcentration.MULT, y = rate, colour = mutant.effect, group = mutant.effect)) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.min, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = rrr.data, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(growth.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
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

IC <- by(growth.data, list(growth.data$mutant.effect, growth.data$cell.type), function(df) {
    sol <- fitSpline(x = df$drug.outsideConcentration.MULT, y = df$rate, pred.y = c(0.1 * rrr.params[[drug]]$psi.max, 0.01 * rrr.params[[drug]]$psi.max, 0), DOF = 15)
    data.frame(
        cell.type = unique(df$cell.type),
        mutant.effect = unique(df$mutant.effect),
        psi90 = sol[1],
        psi99 = sol[2],
        psiZMIC = sol[3]
    )
})
IC <- do.call(rbind, IC)
