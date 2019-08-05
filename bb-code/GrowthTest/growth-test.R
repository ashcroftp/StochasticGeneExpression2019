#' growth-test.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Data directory and cell types to consider
data.directory <- "GrowthTest/zz-data/"
cell.types <- c("REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
drug <- "ciprofloxacin"

## Read summed data (or generate the file first if missing), and plot time-kill curves ----
summed.file <- paste0(data.directory, "summed.dat")
if (!file.exists(summed.file) ) {
    summed <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(cell.type, "growthTest", sep = "_")
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
summed.data$drug <- log10(as.numeric(levels(summed.data$drug.outsideConcentration.MULT))[summed.data$drug.outsideConcentration.MULT])

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

ggplot(summed.data, aes(x = sample.time, y = num.cells, group = drug.outsideConcentration.MULT, colour = drug)) +
    
    geom_hline(data = summed.data[summed.data$sample.time == 0, ], aes(yintercept = num.cells), colour = palette_OkabeIto[8]) +
    
    #geom_vline(data = summed.data[summed.data$sample.time %in% seq(0, max(summed.data$sample.time), 60), ], aes(xintercept = sample.time), colour = palette_OkabeIto[8]) +
    geom_line(data = max.growth, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = max.growth.90, linetype = "dashed", colour = palette_OkabeIto[8]) +
    
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_colour_viridis_c(limits = c(0, 1), name = "drug concentration (log10)") +
    #geom_line(aes(y = fit, colour = drug.dose)) +
    #scale_color_manual(values = gradientPalette(length(levels(summed.data$drug.outsideConcentration.MULT)), palette = "Blues"), name = "drug concentration (x MIC)") +
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


## Load (or first compute) growth rates, and plot these curves ----
growth.file <- paste0(data.directory, "growth-rates.dat")
if (!file.exists(growth.file) ) {
    growth.data <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(cell.type, "growthTest", sep = "_")
        growth <- computeGrowth(data.directory, file.id, final.times = seq(30,300,10))
        growth$cell.type <- factor(cell.type, levels = cell.types)
        return(growth)
    })
    #' Export
    write.csv(growth.data, growth.file, row.names = FALSE)
}

growth.data <- read.csv(growth.file)
growth.data$cell.type <- factor(growth.data$cell.type, levels = cell.types)
growth.data$t.end <- factor(growth.data$t.end, levels = sort(unique(growth.data$t.end)))

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 0.75)
)
rrr.data <- with(
    c(rrr.params[[drug]], list(dose = 10^seq(log10(min(growth.data$drug.outsideConcentration.MULT)), log10(max(growth.data$drug.outsideConcentration.MULT)), length.out = 100))),
    data.frame(cell.type = "WT", mutant.effect = 1, drug.outsideConcentration.MULT = dose, rate =  psi.max - ((psi.max - psi.min) * dose^kappa)/(dose^kappa - psi.min/psi.max), t.end = factor(60) )
)


ggplot(growth.data, aes(x = drug.outsideConcentration.MULT, y = rate, colour = t.end, group = t.end)) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = 0.1 * rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.min, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = rrr.data, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(growth.data$t.end)), palette = "Blues"), name = "duration (minutes)") +
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

## Compare the ICX values of growth rate with those form lineage survival ----
IC <- by(growth.data, list(growth.data$t.end, growth.data$cell.type), function(df) {
    sol <- fitSpline(x = df$drug.outsideConcentration.MULT, y = df$rate, pred.y = c(0.1 * rrr.params[[drug]]$psi.max, 0.01 * rrr.params[[drug]]$psi.max, 0), DOF = 15)
    data.frame(
        cell.type = unique(df$cell.type),
        t.end = unique(df$t.end),
        psi90 = sol[1],
        psi99 = sol[2],
        psiZMIC = sol[3]
    )
})
IC <- do.call(rbind, IC)
#write.csv(IC, file = paste0(data.directory, "MIC.dat"), row.names = F)


## Now compare outputs ICX values from survival and growth data ----
IC <- list(
    survival = read.csv("PopulationSurvival/zz-data/ciprofloxacin_MIC.dat", colClasses = c("factor", "factor", "numeric", "numeric", "numeric")),
    growth = read.csv("GrowthTest/zz-data/MIC.dat", colClasses = c("factor", "factor", "numeric", "numeric", "numeric"))
)

#' #' Normalise relative to WT
#' for (i in c("IC50", "IC90", "IC99")) {
#'     IC$survival[, i] <- IC$survival[, i] / IC$survival[IC$survival$cell.type == "WT", i]
#' }
#' for (i in c("psi90", "psi99", "psiZMIC")) {
#'     for (j in levels(IC$growth$t.end)) {
#'         IC$growth[IC$growth$t.end == j, i] <- IC$growth[IC$growth$t.end == j, i] / IC$growth[IC$growth$t.end == j & IC$growth$cell.type == "WT", i]
#'     }
#' }
#' rm(i)

#' Melt each data.frame
IC$survival$cell.type <- factor(IC$survival$cell.type, levels = cell.types)
IC$survival <- IC$survival[IC$survival$mutant.effect %in% c("1", "200"), ]
IC$survival$mutant.effect <- factor(IC$survival$mutant.effect, levels =  c("1", "200"))
IC$survival <- melt(IC$survival, id.vars = c("cell.type", "mutant.effect"))

IC$growth$cell.type <- factor(IC$growth$cell.type, levels = cell.types)
IC$growth <- melt(IC$growth, id.vars = c("cell.type", "t.end"))
IC$growth$t.end <- as.numeric(levels(IC$growth$t.end)[IC$growth$t.end])
levels(IC$growth$variable) <- list("10% max" = "psi90", "1% max" = "psi99", "zero growth" = "psiZMIC")


ggplot(IC$growth, aes(x = t.end, y = value, shape = variable, colour = variable)) +
    geom_hline(data = IC$survival, aes(yintercept = value, linetype = variable), colour = "black") +
    geom_point(size = 1.0, stroke = 0.5, fill = "white") +
    facet_wrap(~cell.type) +
    #scale_y_log10(breaks = logTicks()$major, minor_breaks = logTicks()$minor) +
    scale_y_log10(breaks = c(seq(0.7,0.9,0.1), seq(1,10,1)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    scale_colour_manual(values = gradientPalette(length(levels(IC$growth$variable)), palette = "Oranges"), name = "growth rate\nmeasure of IC") +
    scale_shape_manual(values = c(21,23,24), name = "growth rate\nmeasure of IC") +
    scale_linetype_manual(values = c(IC50 = "solid", IC90 = "dashed", IC99 = "dotted"), name = "lineage survival\nmeasure of IC") +
    labs(x = "experiment duration (minutes)", y = "inhibitory concentration (x WT MIC)") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.66,0.15),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")
    )

# IC <- do.call(rbind, list(IC$survival[, c("group", "cell.type", "measure", "value")], IC$growth[, c("group", "cell.type", "measure", "value")]))
# IC <- IC[!IC$cell.type == "WT", ]
# IC$cell.type <- factor(IC$cell.type, levels = c("REG-ON", "REG-OFF", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO"))
# 
# ggplot(IC, aes(x = measure, y = value, colour = cell.type)) +
#     geom_hline(yintercept = 1) +
#     geom_point(position = position_dodge(0.5)) +
#     scale_color_manual(values = palette_OkabeIto[1:7], name = "cell type") +
#     facet_grid(~ group, scales = "free_x", space = "free_x", labeller = label_panels) +
#     labs(x = "measure to determine MIC", y = "MIC relative to WT") +
#     theme_bw() +
#     theme(
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom"
#     )
