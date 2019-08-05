#' plot-WT-molecules.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:
source("R/basicFunctions.R")
source("R/plotFunctions.R")
source("Simulate/filenames.R")
library("plyr")
library("MASS")


#' Load measurement data
cell.type <- "WT"
data.directory <- "MoleculeDistributions/zz-data/"
file.id <- paste(cell.type, "molecules", sep = "_")
molecules.data <- read.csv(output.file(data.directory, file.id))

#' Get fitted distribution
efflux.range <- c(0:max(molecules.data$efflux.protein))
fit <- fitdistr(x = molecules.data$efflux.protein, "negative binomial")$estimate
fit["prob"] <- fit[["mu"]]/(fit[["size"]] + fit[["mu"]])
efflux.dist <- data.frame(number = efflux.range, density = dnbinom(efflux.range, size = fit[["size"]], mu = fit[["mu"]]))
#' Parameters of the corresponding gamma distribution
gamma.pars <- list(shape = fit[["size"]]*fit[["prob"]], rate = 1 - fit[["prob"]], burst.size = 1/(1 - fit[["prob"]]))

#' Plot
density.plot <- function(data, variable, colour) {
    data <- data.frame(number = data[, variable])
    if (variable == "efflux.DNA") {
        data$number <- factor(data$number, levels = c("0", "1"))
        levels(data$number) <- c(`0` = "off", `1` = "on")
    }
    counts <- count(data)
    names(counts) <- c("number", "count")
    counts$density <- counts$count/sum(counts$count)

    if (variable == "efflux.DNA") {
        plot <- ggplot(counts, aes(x = number, y = density)) +
            geom_col(fill = colour) +
            theme_bw() +
            theme(
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank()
            ) +
            labs(x = "DNA state", y = "density")
    }
    else {
        plot <- ggplot(counts, aes(x = number, y = density)) +
            #geom_segment(aes(xend = number), yend = -4, lineend = "square", colour = colour) +
            geom_point(colour = colour) +
            geom_vline(xintercept = mean(data$number), linetype = "dashed") +
            #geom_bar(stat = "identity") +
            theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            labs(x = ifelse(variable == "efflux.mRNA", "number of mRNA", ifelse(variable == "efflux.protein", "number of efflux proteins", "number of target proteins")), y = "density")
    }
    return(plot)
}

plot.list <- list(
    efflux.DNA = density.plot(molecules.data, "efflux.DNA", colour = palette_OkabeIto[2]),
    efflux.mRNA = density.plot(molecules.data, "efflux.mRNA", colour = palette_OkabeIto[2]) +
        scale_y_log10(labels = sciFormat),
    efflux.protein = density.plot(molecules.data, "efflux.protein", colour = palette_OkabeIto[2]) +
        geom_line(data = efflux.dist, colour = palette_OkabeIto[8]) +
        scale_y_log10(labels = sciFormat, limits = c(1e-4, 1e0)),
    target.protein = density.plot(molecules.data, "target.protein", colour = palette_OkabeIto[1])
)
plot_grid(plotlist = plot.list, axis = "lb", labels = c("AUTO"), label_size = 12, label_fontface = "plain")
