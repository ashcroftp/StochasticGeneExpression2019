#' systematic.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/filenames.R")

#' Data directory and cell types to consider
data.directory <- "Systematic/zz-data/"
cell.types <- c("WT", "REG-ON", "STRUCT-BIND", "STRUCT-CAT")

## 1. Compute survival probabilities ----
#' If survival probability file does not exist, we need to compute this form the raw simulation data
survival.file <- paste0(data.directory, "survival.dat")
if (!file.exists(survival.file) ) {
    survival <- DFapply(cell.types, function(cell.type) {
        file.id <- paste("systematic", cell.type, sep = "_")
        survival.df <- computeSurvival(data.directory, file.id)
        #' Add cell type to DF
        survival.df$cell.type <- factor(cell.type, levels = cell.types)
        return(survival.df)
    })
    #' Export
    write.csv(survival, survival.file, row.names = FALSE)
}


#' If the simplified survival file doesn't exist, then we need to compute it too
simple.survival.file <- paste0(data.directory, "simple-survival.dat")
if (!file.exists(simple.survival.file) ) {
    #' Load the survival file
    survival <- read.csv(survival.file)
    #' Now look up which cell type (1-729) each row of these results corresponds too,
    #' so that we can create a cleaner data frame for comparisons
    params.output <- read.csv(paste0(data.directory, "params.output"))
    params.output$id <- factor(seq_len(nrow(params.output)))
    #' Truncate numerical outputs for efficient matching
    for (param in names(survival)[names(survival) %in% names(params.output)]) {
        params.output[, param] <- signif(params.output[, param], 7)
        survival[, param] <- signif(survival[, param], 7)
    }
    rm(param)
    
    survival.simple <- merge(params.output, survival, all.y = TRUE)
    survival.simple <- survival.simple[, c("cell.type", "id", "binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC", "MIC.fraction", "mutant.effect", "drug.outsideConcentration.MULT", "survive", "extinction.time", "extinction.time.sd", "num.reps")]
    #' Export
    write.csv(survival.simple, simple.survival.file, row.names = FALSE)
}

## 2. Plot the survival curves of a given cell type (or a list of them) ----
#' Load the data
data <- read.csv(simple.survival.file)
#' Refactor
data <- data[data$cell.type %in% cell.types, ]
data$cell.type <- factor(data$cell.type, levels = cell.types)
data$id <- factor(data$id, levels = sort(unique(data$id)))
#data$param.index <- factor(data$param.index, levels = sort(unique(data$param.index)))
data$mutant.effect <- factor(data$mutant.effect, levels = sort(unique(data$mutant.effect)))

#' Plot out IC spline fit for a given cell type
plotKillCurve <- function(data, cell.type, id, S.levels = c(0.5, 0.1)) {
    df <- data[data$cell.type == cell.type & data$id == id, c("drug.outsideConcentration.MULT", "survive", "num.reps")]
    names(df) <- c("x", "y", "weight")

    data.plot <- ggplot(df, aes(x = x, y = y)) +
        geom_hline(yintercept = S.levels) +
        #scale_y_log10() +
        geom_point()

    my.spline <- fitSpline(x = df$x, y = df$y, weights = df$weight, DOF = 15)
    IC.points <- data.frame(x = fitSpline(df$x, y = df$y, weights = df$weight, pred.y = S.levels), y = S.levels, DOF = 15)
    IC.points.lin <- data.frame(x = fitLinear(df$x, y = df$y, pred.y = S.levels), y = S.levels)
    #my.spline <- fitSigmoid(x = df$x, y = df$y)
    #IC.points <- data.frame(x = fitSigmoid(df$x, y = df$y, pred.y = c(0.5,0.1,0.01,0.001)), y = c(0.5,0.1,0.01,0.001))

    data.plot + geom_line(data = my.spline) + geom_point(data = IC.points, colour = "red") + geom_point(data = IC.points.lin, colour = "pink", shape = 4)
}
#plotKillCurve(data, "STRUCT-BIND", "727", S.levels = c(0.5,0.1,0.01))
plot_grid(plotlist = lapply(c(108,169,178,196,197,198,206,207,215,216,439,466,467,701), function(id) plotKillCurve(data, "STRUCT-BIND", as.character(id), S.levels = c(0.5,0.1,0.01))))
plot_grid(plotlist = lapply(c(163,164,165,172,173,174,181,182,183,678,687,696), function(id) plotKillCurve(data, "STRUCT-CAT", as.character(id), S.levels = c(0.5,0.1,0.01))))


## 3. Compute the ICX values and save the data ----
IC.file <- paste0(data.directory, "IC.dat")
if (!file.exists(IC.file) ) {
    #' Load the data
    data <- read.csv(simple.survival.file)
    #' Refactor
    data <- data[data$cell.type %in% cell.types, ]
    data$cell.type <- factor(data$cell.type, levels = cell.types)
    data$id <- factor(data$id, levels = sort(unique(data$id)))
    #data$param.index <- factor(data$param.index, levels = sort(unique(data$param.index)))
    data$mutant.effect <- factor(data$mutant.effect, levels = sort(unique(data$mutant.effect)))
    
    IC <- by(data, list(data$id, data$cell.type), function(df) {
        #' Order df by drug concentration
        df <- df[order(df$drug.outsideConcentration.MULT), ]
        #' Extract the ICx values
        #sol <- computeICX(df)
        sol <- fitSpline(df$drug.outsideConcentration.MULT, y = df$survive, weights = df$num.reps, pred.y = 1 - c(0.5, 0.9, 0.99), DOF = 15)
        
        #' else
        out <- data.frame(
            cell.type = unique(df$cell.type),
            #mutant.effect = unique(df$mutant.effect),
            id = unique(df$id),
            IC50 = sol[1],
            IC90 = sol[2],
            IC99 = sol[3]
        )
    })
    IC <- do.call(rbind, IC)
    IC[is.na(IC$IC99), ]
    paste(IC[IC$cell.type == "STRUCT-BIND" & is.na(IC$IC99), "id"], collapse = ",")
    #' Export
    write.csv(IC, IC.file, row.names = FALSE)
}

## 4. Plots of the ICX values across the different cell types ----
IC <- read.csv(paste0(data.directory, "IC.dat"))
IC$cell.type <- factor(IC$cell.type, levels = cell.types)
IC$id <- factor(IC$id, levels = sort(unique(IC$id)))

IC.melt <- melt(IC, id.vars = c("cell.type", "id"))

ggplot(IC.melt, aes(x = cell.type, y = value, group = id)) +
    geom_point(alpha = 0.2) +
    geom_line(alpha = 0.2) +
    scale_y_log10() +
    facet_wrap(~variable, labeller = label_panels) +
    labs(x = NULL, y = "inhibitory concentration") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

ggplot(IC.melt, aes(x = cell.type, y = value, fill = cell.type, colour = cell.type)) +
    geom_hline(yintercept = 1, colour = palette_OkabeIto[8]) +
    geom_violin() +
    scale_fill_manual(values = cellColours(), name = "mutation") +
    scale_colour_manual(values = cellColours(), name = "mutation") +
    scale_y_log10() +
    facet_wrap(~variable, labeller = label_panels) +
    labs(x = NULL, y = "inhibitory concentration") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    )



## 5. Now we recombine the IC dataframe with the parameters, and look for correlations between parameters and the measured ICX values ----
params.output <- read.csv(paste0(data.directory, "params.output"))
params.output$id <- factor(seq_len(nrow(params.output)))
params.output <- params.output[, c("id", "binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC", "MIC.fraction")]

IC <- read.csv(paste0(data.directory, "IC.dat"))
IC$cell.type <- factor(IC$cell.type, levels = cell.types)
IC$id <- factor(IC$id, levels = sort(unique(IC$id)))

IC <- merge(IC, params.output, by = "id")

#' Inspect the results vs variables
#' First we melt the variables
melted <- melt(
    data = IC,
    id.vars = c("id", "cell.type", "IC50", "IC90", "IC99"),
    measure.vars = c("binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC", "MIC.fraction"),
    variable.name = "variable", value.name = "variable.value"
)
#' Now we melt the response variables
melted <- melt(
    data = melted,
    id.vars = c("id", "cell.type", "variable", "variable.value"),
    measure.vars = c("IC50", "IC90", "IC99"),
    variable.name = "IC", value.name = "IC.value"
)

ggplot(melted, aes(x = variable.value, y = IC.value)) +
    geom_smooth(method = "lm", se = F, colour = palette_OkabeIto[8]) +
    geom_point() +
    #scale_color_viridis_c(name = "MIC fraction", limits = c(0,1)) +#(trans = "log10") +
    facet_grid(IC + cell.type ~ variable, scales = "free_x", margins = c("cell.type", "IC")) +
    scale_x_log10() +
    scale_y_log10() + 
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )


#' Compute changes in IC relative to WT
delta.IC <- by(IC, list(IC$id), function(df) {
    for (i in c("IC50", "IC90", "IC99")) {
        df[, i] <- df[, i] / df[df$cell.type == "WT", i]
    }
    return(df)
})
delta.IC <- do.call(rbind, delta.IC)
delta.IC.melt <- melt(delta.IC, id.vars = c("cell.type", "id", "binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC", "MIC.fraction"))

ggplot(delta.IC.melt, aes(x = cell.type, y = value, group = id)) +
    geom_point(alpha = 0.2) +
    geom_line(alpha = 0.2) +
    scale_y_log10() +
    facet_wrap(~variable, labeller = label_panels) +
    labs(x = NULL, y = "inhibitory concentration") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

ggplot(delta.IC.melt, aes(x = cell.type, y = value, fill = cell.type, colour = cell.type)) +
    geom_hline(yintercept = 1, colour = palette_OkabeIto[8]) +
    geom_violin() +
    scale_fill_manual(values = cellColours(), name = "mutation") +
    scale_colour_manual(values = cellColours(), name = "mutation") +
    scale_y_log10() +
    facet_wrap(~variable, labeller = label_panels) +
    labs(x = NULL, y = "inhibitory concentration") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    )


ggplot(delta.IC.melt[delta.IC.melt$efflux.pumps > delta.IC.melt$targets, ], aes(x = cell.type, y = value, fill = cell.type, colour = cell.type)) +
    geom_hline(yintercept = 1, colour = palette_OkabeIto[8]) +
    geom_violin() +
    scale_fill_manual(values = cellColours(), name = "mutation") +
    scale_colour_manual(values = cellColours(), name = "mutation") +
    scale_y_log10() +
    facet_wrap(~variable, labeller = label_panels) +
    labs(x = NULL, y = "inhibitory concentration") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top"
    )

#' Inspect the results vs variables
#' First we melt the variables
delta.melted <- melt(
    data = delta.IC,
    id.vars = c("id", "cell.type", "IC50", "IC90", "IC99"),
    measure.vars = c("binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC", "MIC.fraction"),
    variable.name = "variable", value.name = "variable.value"
)
#' Now we melt the response variables
delta.melted <- melt(
    data = delta.melted,
    id.vars = c("id", "cell.type", "variable", "variable.value"),
    measure.vars = c("IC50", "IC90", "IC99"),
    variable.name = "IC", value.name = "IC.value"
)


plot_grid(
    plotlist = lapply(cell.types[cell.types != "WT"], function(cell.type) {
        ggplot(delta.melted[delta.melted$cell.type == cell.type, ], aes(x = variable.value, y = IC.value)) +
            geom_smooth(method = "lm", se = F, colour = palette_OkabeIto[8]) +
            geom_point() +
            #scale_color_viridis_c(name = "MIC fraction", limits = c(0,1)) +#(trans = "log10") +
            facet_grid(IC ~ variable, scales = "free_x") +
            scale_x_log10() + #labels = sciFormat) +
            scale_y_log10() + #labels = sciFormat) + 
            theme_bw() +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(hjust = 0),
                #panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "None"
            )
    }),
    labels = cell.types[cell.types != "WT"]
)

ggplot(delta.melted[delta.melted$cell.type != "WT", ], aes(x = variable.value, y = IC.value, colour = cell.type)) +
    geom_smooth(method = "lm", se = F) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = cellColours()) +
    #scale_color_viridis_c(name = "MIC fraction", limits = c(0,1)) +#(trans = "log10") +
    facet_grid(IC ~ variable, scales = "free_x") +
    scale_x_log10() + #labels = sciFormat) +
    scale_y_log10() + #labels = sciFormat) + 
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

facet.labels <- c(binding = "binding rate", catalysis = "catalysis rate", diffusion = "diffusion rate", targets = "ave. num. targets per cell", efflux.pumps = "ave. num. efflux pumps per cell", drug.molecules = "num. drug molecules per cell @ MIC", MIC = "external drug conc. @ MIC", MIC.fraction = "fraction of bound targets @ MIC")
ggplot(delta.melted[delta.melted$cell.type != "WT" & delta.melted$IC == "IC90", ], aes(x = variable.value, y = IC.value, colour = cell.type)) +
    geom_smooth(method = "lm", se = F) +
    geom_point(alpha = 0.5) +
    scale_colour_manual(values = cellColours(), name = "mutation") +
    #scale_color_viridis_c(name = "MIC fraction", limits = c(0,1)) +#(trans = "log10") +
    facet_grid(~ variable, scales = "free_x", switch = "x", labeller = as_labeller(facet.labels)) +
    scale_x_log10(labels = sciFormat) +
    scale_y_log10() + #labels = sciFormat) +
    labs(x = NULL, y = "IC90 (relative to WT)") +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

getPlot <- function(data) {
    library(gtable)
    library(grid)
    
    facet.labels <- c(binding = "binding rate", catalysis = "catalysis rate", diffusion = "diffusion rate", targets = "ave. num. targets per cell", efflux.pumps = "ave. num. efflux pumps per cell", drug.molecules = "num. drug molecules per cell @ MIC", MIC = "external drug conc. @ MIC", MIC.fraction = "fraction of bound targets @ MIC")
    p <- ggplot(data, aes(x = variable.value, y = IC.value, colour = cell.type)) +
        geom_smooth(method = "lm", se = F) +
        geom_point(alpha = 0.5, position = position_dodge(width = 0.2)) +
        #geom_jitter(alpha = 0.5, height = 0, width = 0.2) +
        scale_colour_manual(values = cellColours(), name = "mutation") +
        #scale_color_viridis_c(name = "MIC fraction", limits = c(0,1)) +#(trans = "log10") +
        facet_grid(~ variable, scales = "free_x", labeller = as_labeller(facet.labels)) +
        scale_x_log10(labels = sciFormat) +
        scale_y_log10() + #labels = sciFormat) +
        labs(x = NULL, y = "IC90 (relative to WT)") +
        theme_bw() +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0.5),
            #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom"
        )
    
    # Convert the plot to a grob
    gt <- ggplotGrob(p)
    
    # Get the positions of the panels in the layout: t = top, l = left, ...
    panels <- c(subset(gt$layout, grepl("panel", gt$layout$name), select = t:r))
    
    # Add a row below the x-axis tick mark labels,
    # the same height as the strip
    gt = gtable_add_rows(gt, gt$height[min(panels$t) - 1], max(panels$b) + 2)
    
    # Get the strip grob
    stripGrob = gtable_filter(gt, "strip-t")
    
    # Insert the strip grob into the new row
    gt = gtable_add_grob(gt, stripGrob, t = max(panels$b) + 3, l = min(panels$l), r = max(panels$r))
    
    # remove the old strip
    gt = gt[-(min(panels$t) - 1), ]
    
    grid.newpage()
    grid.draw(gt)
}
getPlot(delta.melted[delta.melted$cell.type != "WT" & delta.melted$IC == "IC90", ])

getPlot(delta.melted[delta.melted$cell.type %in% c("STRUCT-BIND", "STRUCT-CAT") & delta.melted$IC == "IC90", ])
getPlot(delta.melted[delta.melted$cell.type %in% c("STRUCT-BIND", "STRUCT-CAT") & delta.melted$IC == "IC90" & delta.melted$id %in% params.output[params.output$efflux.pumps > params.output$targets, "id"], ])


## 6. Correlations between different mutants ----
params.output <- read.csv(paste0(data.directory, "params.output"))
params.output$id <- factor(seq_len(nrow(params.output)))
params.output <- params.output[, c("id", "binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC", "MIC.fraction")]

IC <- read.csv(paste0(data.directory, "IC.dat"))
IC$cell.type <- factor(IC$cell.type, levels = cell.types)
IC$id <- factor(IC$id, levels = sort(unique(IC$id)))

IC <- dcast(IC, id ~ cell.type, value.var = "IC90")

for(i in cell.types) IC[,i] <- log10(IC[,i])

library(GGally)
ggpairs(IC, columns = cell.types, )
