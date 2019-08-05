#' efflux-bias-test.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Cell types and data directory
cell.types <- c("WT", "REG-ON", "STRUCT-BIND")
data.directory <- "EffluxBias/zz-data/"
drug <- "ciprofloxacin"

## Survival probability ----

#' If the output file doesn't exist, we need to compute the survival probability
survival.file <- paste0(data.directory, "survival.dat")
if (!file.exists(survival.file) ) {
    survival.df <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(cell.type, "survival", sep = "_")
        survival <- computeSurvival(data.directory, file.id)
        survival$cell.type <- factor(cell.type, levels = cell.types)
        return(survival)
    })
    #' Export
    write.csv(survival.df, survival.file, row.names = FALSE)
}

#' Load the data
data <- read.csv(survival.file)
#' Refactor
data <- data[data$cell.type %in% cell.types, ]
data$cell.type <- factor(data$cell.type, levels = cell.types)
data$param.index <- factor(data$param.index, levels = sort(unique(data$param.index)))
data$mutant.effect <- factor(data$mutant.effect, levels = sort(unique(data$mutant.effect)))
data$efflux.binom.param <- factor(data$efflux.binom.param, levels = sort(unique(data$efflux.binom.param)))
#' Plot the survival probabilities on a log-x scale
# #survival.plot <- 
# ggplot(data[data$mutant.effect %in% c("1","200"), ], aes(x = drug.outsideConcentration.MULT, y = survive, colour = efflux.binom.param, group = efflux.binom.param)) +
#     geom_vline(xintercept = 1, colour = palette_OkabeIto[8]) +
#     geom_line(size = 1.2) +
#     geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
#     scale_color_manual(values = colorRampPalette(brewer.pal(11, "PiYG")[5:1])(length(levels(data$efflux.binom.param))), name = "efflux bias") +
#     #scale_color_manual(values = gradientPalette(length(levels(data$efflux.binom.param)), palette = "Purples"), name = "efflux bias") +
#     scale_x_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,9,1), seq(10,100,10)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
#     facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
#     theme_bw() +
#     theme(
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         legend.position = "top"
#     ) +
#     labs(x = "drug concentration (x MIC)", y = "survival probability")
# 
# data$drug <- factor(data$drug.outsideConcentration.MULT)
# ggplot(data, aes(x = drug, y = efflux.binom.param)) +
#     geom_tile(aes(fill = survive)) +
#     scale_fill_distiller(palette = "Purples", name = "survival probability", type = 'seq', na.value = "grey60", direction = 1, limits = c(0,1)) +
#     #scale_fill_viridis_c(name = "survival probability", direction = 1, limits = c(0,1)) +
#     facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
#     coord_fixed() +
#     theme_bw() +
#     theme(
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0)
#     ) +
#     labs(x = "drug concentration (x MIC)", y = "efflux bias") +
#     scale_y_discrete(expand = c(0,0)) +
#     scale_x_discrete(expand = c(0,0))
# 
# ggplot(data[data$mutant.effect %in% c("1","200"), ], aes(x = as.numeric(levels(efflux.binom.param)[efflux.binom.param]), y = survive, colour = drug, group = drug)) +
#     geom_vline(xintercept = 1, colour = palette_OkabeIto[8]) +
#     geom_line(size = 1.2) +
#     geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
#     scale_color_manual(values = colorRampPalette(brewer.pal(11, "PiYG")[5:1])(length(levels(data$drug))), name = "drug concentration (x MIC)") +
#     #scale_color_manual(values = gradientPalette(length(levels(data$efflux.binom.param)), palette = "Purples"), name = "efflux bias") +
#     facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
#     theme_bw() +
#     theme(
#         strip.background = element_blank(),
#         strip.text = element_text(hjust = 0),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         legend.position = "top"
#     ) +
#     labs(x = "efflux bias", y = "survival probability")

#' Compute IC90 and plot this
IC <- by(data, list(data$efflux.binom.param, data$mutant.effect, data$cell.type), function(df) {
    sol <- fitSpline(x = df$drug.outsideConcentration.MULT, y = df$survive, weights = df$num.reps, pred.y = c(0.5, 0.1, 0.01, 0.001), DOF = 15)
    data.frame(
        cell.type = unique(df$cell.type),
        mutant.effect = unique(df$mutant.effect),
        efflux.bias = unique(df$efflux.binom.param),
        IC50 = sol[1],
        IC90 = sol[2],
        IC99 = sol[3],
        IC99.9 = sol[4]
    )
})
IC <- do.call(rbind, IC)

#survival.plot <-
ggplot(IC, aes(x = as.numeric(levels(efflux.bias)[efflux.bias]), y = IC90, colour = mutant.effect, group = mutant.effect)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    scale_y_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,9,1), seq(10,100,10)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    facet_wrap(~cell.type, ncol = 3, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "efflux bias", y = "inhibitory concentration (IC90)")

## Read summed data (or generate the file first if missing), and plot time-kill curves ----
summed.file <- paste0(data.directory, "summed.dat")
if (!file.exists(summed.file) ) {
    summed <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(cell.type, "growth", sep = "_")
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
summed.data$mutant.effect <- factor(summed.data$mutant.effect, levels = sort(unique(summed.data$mutant.effect)))
summed.data$efflux.binom.param <- factor(summed.data$efflux.binom.param, levels = sort(unique(summed.data$efflux.binom.param)))

## Compute growth rates from the summed data and plot dose-response curves ----
growth.file <- paste0(data.directory, "growth.dat")
if (!file.exists(growth.file) ) {
    growth.df <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(cell.type, "growth", sep = "_")
        growth.rates <- computeGrowth(data.directory, file.id, final.times = 180)
        growth.rates$cell.type <- factor(cell.type, levels = cell.types)
        return(growth.rates)
    })
    #' Export
    write.csv(growth.df, growth.file, row.names = FALSE)
}
growth.data <- read.csv(growth.file)
growth.data$cell.type <- factor(growth.data$cell.type, levels = cell.types)
growth.data$mutant.effect <- factor(growth.data$mutant.effect, levels = sort(unique(growth.data$mutant.effect)))
growth.data$efflux.binom.param <- factor(growth.data$efflux.binom.param, levels = sort(unique(growth.data$efflux.binom.param)))

#' Parameters from Regoes et al. (2004) for comparison of growth rates
rrr.params <- list(
    ciprofloxacin = list(psi.max = 0.88/60, psi.min = -6.5/60, kappa = 1.1),
    rifampicin = list(psi.max = 0.70/60, psi.min = -4.3/60, kappa = 0.75)
)
rrr.data <- with(
    c(rrr.params[[drug]], list(dose = 10^seq(log10(min(growth.data$drug.outsideConcentration.MULT)), log10(max(growth.data$drug.outsideConcentration.MULT)), length.out = 100))),
    data.frame(cell.type = "WT", mutant.effect = 1, efflux.binom.param = factor(0.5), drug.outsideConcentration.MULT = dose, rate =  psi.max - ((psi.max - psi.min) * dose^kappa)/(dose^kappa - psi.min/psi.max) )
)

ggplot(growth.data[growth.data$mutant.effect %in% c("1","200"), ], aes(x = drug.outsideConcentration.MULT, y = rate, colour = efflux.binom.param, group = efflux.binom.param)) +
    geom_hline(yintercept = 0, colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.max, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_hline(yintercept = rrr.params[[drug]]$psi.min, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(data = rrr.data, linetype = "dashed", colour = palette_OkabeIto[8]) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    #scale_color_manual(values = gradientPalette(length(levels(growth.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    scale_color_manual(values = colorRampPalette(brewer.pal(11, "PiYG")[5:1])(length(levels(growth.data$efflux.binom.param))), name = "efflux bias") +
    
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

## Compute ICX values ----
zMIC <- by(growth.data, list(growth.data$efflux.binom.param, growth.data$mutant.effect, growth.data$cell.type), function(df) {
    sol <- fitSpline(x = df$drug.outsideConcentration.MULT, y = df$rate, pred.y = c(0.1 * rrr.params[[drug]]$psi.max, 0.01 * rrr.params[[drug]]$psi.max, 0), DOF = 15)
    data.frame(
        cell.type = unique(df$cell.type),
        mutant.effect = unique(df$mutant.effect),
        efflux.bias = unique(df$efflux.binom.param),
        psi90 = sol[1],
        psi99 = sol[2],
        psiZMIC = sol[3]
    )
})
zMIC <- do.call(rbind, zMIC)

ggplot(zMIC, aes(x = as.numeric(levels(efflux.bias)[efflux.bias]), y = psiZMIC, colour = mutant.effect, group = mutant.effect)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(zMIC$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    scale_y_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,9,1), seq(10,100,10)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "efflux bias", y = "inhibitory concentration (zMIC)")


## Combine
dat <- rbind(
    data.frame(cell.type = IC$cell.type, mutant.effect = IC$mutant.effect, efflux.bias = IC$efflux.bias, ic = IC$IC90, experiment = factor("IC90")),
    data.frame(cell.type = zMIC$cell.type, mutant.effect = zMIC$mutant.effect, efflux.bias = zMIC$efflux.bias, ic = zMIC$psiZMIC, experiment = factor("zMIC"))
)

ggplot(dat, aes(x = as.numeric(levels(efflux.bias)[efflux.bias]), y = ic, colour = mutant.effect, linetype = experiment)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    scale_y_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,9,1), seq(10,100,10)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    facet_wrap(~cell.type, ncol = 3, labeller = label_panels) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "efflux bias", y = "inhibitory concentration")
