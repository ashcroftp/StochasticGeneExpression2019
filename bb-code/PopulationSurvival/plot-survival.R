#' plot-survival.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Drug and cell types to consider
drug <- "ciprofloxacin"
cell.types <- c("REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT", "WT", "KO")
#drug <- "cefotaxime"
#cell.types <- c("WT", "STRUCT-BIND", "STRUCT-CAT")
data.directory <- "PopulationSurvival/zz-data/"

## Survival probability vs drug concentration ----

#' If the output file doesn't exist, we need to compute the survival probability
survival.file <- paste0(data.directory, paste(drug, "survival.dat", sep = "_"))
if (!file.exists(survival.file) ) {
    survival.df <- DFapply(cell.types, function(cell.type) {
        file.id <- paste(drug, cell.type, sep = "_")
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
data$drug.off <- factor(data$drug.off, levels = sort(unique(data$drug.off)))


## Survival probabilities and extinction times ----
survival.data <- data[data$drug.off == Inf, ]

#' #' Plot the survival probabilities
#' survival.plot <- ggplot(survival.data, aes(x = drug.outsideConcentration.MULT, y = survive, group = mutant.effect)) +
#'     geom_vline(xintercept = 1, colour = palette_OkabeIto[8]) +
#'     geom_line(aes(colour = mutant.effect), size = 1.2) +
#'     geom_point(aes(colour = mutant.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
#'     scale_color_manual(values = gradientPalette(length(levels(survival.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
#'     facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
#'     theme_bw() +
#'     theme(
#'         strip.background = element_blank(),
#'         strip.text = element_text(hjust = 0),
#'         panel.grid.minor = element_blank(),
#'         legend.position = "top"
#'     ) +
#'     labs(x = "drug concentration (x MIC)", y = "survival probability")

#' Plot the survival probabilities on a log-x scale
survival.plot <- ggplot(survival.data, aes(x = drug.outsideConcentration.MULT, y = survive, group = mutant.effect)) +
    geom_vline(xintercept = 1, colour = palette_OkabeIto[8]) +
    geom_line(aes(colour = mutant.effect), size = 1.2) +
    geom_point(aes(colour = mutant.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    scale_x_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,10,1)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    coord_cartesian(xlim = c(0.4,10)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "top"
    ) +
    labs(x = "drug concentration (x MIC)", y = "survival probability")

## Compute the ICX values from this survival data ----
IC <- by(survival.data, list(survival.data$mutant.effect, survival.data$cell.type), function(df) {
    sol <- fitSpline(x = df$drug.outsideConcentration.MULT, y = df$survive, weights = df$num.reps, pred.y = c(0.5, 0.1, 0.01), DOF = 20)
    data.frame(
        cell.type = unique(df$cell.type),
        mutant.effect = unique(df$mutant.effect),
        IC50 = sol[1],
        IC90 = sol[2],
        IC99 = sol[3]
    )
})
IC <- do.call(rbind, IC)
#write.csv(IC, file = "PopulationSurvival/zz-data/ciprofloxacin_MIC.dat", row.names = F)

# spline.pred <- by(survival.data, list(survival.data$mutant.effect, survival.data$cell.type), function(df) {
#     sol <- fitSpline(x = df$drug.outsideConcentration.MULT, y = df$survive, DOF = 20)
#     data.frame(
#         cell.type = unique(df$cell.type),
#         mutant.effect = unique(df$mutant.effect),
#         drug.outsideConcentration.MULT = sol$x,
#         survive = sol$y
#     )
# })
# spline.pred <- do.call(rbind, spline.pred)
# survival.plot + geom_line(data = spline.pred) + geom_hline(yintercept = c(0.5,0.1,0.01), colour = palette_OkabeIto[8])


#' Plot the extinction times (with standard deviation)
ggplot(survival.data, aes(x = drug.outsideConcentration.MULT, y = extinction.time, group = mutant.effect)) +
    geom_vline(xintercept = 1, colour = palette_OkabeIto[8]) +
    geom_line(aes(colour = mutant.effect), size = 1.2) +
    geom_errorbar(aes(colour = mutant.effect, ymin = extinction.time - extinction.time.sd, ymax = extinction.time + extinction.time.sd), size = 0.1) +
    geom_point(aes(colour = mutant.effect), size = 1.0, shape = 21, stroke = 0.5, fill = "white", show.legend = F) +
    scale_color_manual(values = gradientPalette(length(levels(survival.data$mutant.effect)), palette = "Blues"), name = "mutant effect") +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    scale_x_log10(breaks = c(seq(0.1,0.9,0.1), seq(1,10,1)), labels = function(x) ifelse(x == 0.5, "0.5", ifelse(x == 1.0, "1", ifelse(x > 4.5 & x < 5.5, "5", ifelse(x == 10.0, "10", ""))))) +
    coord_cartesian(xlim = c(0.4,10), ylim = c(0,500)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
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
    scale_fill_distiller(palette = "Purples", name = "survival probability", type = 'seq', na.value = "grey60", direction = 1, limits = c(0,1)) +
    #scale_fill_viridis_c(name = "survival probability", direction = 1, limits = c(0,1)) +
    facet_wrap(~cell.type, ncol = 2, labeller = label_panels) +
    coord_fixed() +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0)
    ) +
    labs(x = "drug concentration (x MIC)", y = "treatment duration (minutes)") +
    scale_y_discrete(expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0))
