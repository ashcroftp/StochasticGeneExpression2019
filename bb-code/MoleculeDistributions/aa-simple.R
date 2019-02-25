#' aa-simple.R
#' Author: Peter Ashcroft, ETH Zurich

#' A simple script to plot the distributions of molecules in cells in the absence of antibiotics.

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/plotFunctions.R")
library(plyr) # For the `count` function

# Simulate ----
#' Load parameters & modify (and make sure drug is set to zero)
basic.params <- loadParameters(drug = "ciprofloxacin")
basic.params$drug.outsideConcentration <- 0

#' Define some initial cell types and then create an explicit vector
num.cells <- c(`KO` = 0, `WT` = 500, `REG-ON` = 0, `REG-OFF` = 0, `REG-BURST` = 0, `STRUCT-BIND` = 0, `STRUCT-CAT` = 0)
#' Set the seed for reproducibility
set.seed(1111)
#' Simulate multiple cells for multiple generstion, and just retain the last timepoint
#' Here we can use the initial-condition function, which does just this
cell.data <- simulateIC(basic.params, num.cells, t.max = 200)

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
            geom_point(colour = colour) +
            geom_vline(xintercept = mean(data$number), linetype = "dashed") +
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
    efflux.DNA = density.plot(cell.data, "efflux.DNA", colour = palette_OkabeIto[2]),
    efflux.mRNA = density.plot(cell.data, "efflux.mRNA", colour = palette_OkabeIto[2]) + scale_y_log10(labels = sciFormat),
    efflux.protein = density.plot(cell.data, "efflux.protein", colour = palette_OkabeIto[2]) + scale_y_log10(labels = sciFormat, limits = c(1e-4, 1e0)),
    target.protein = density.plot(cell.data, "target.protein", colour = palette_OkabeIto[1])
)
plot_grid(plotlist = plot.list, axis = "lb", labels = c("AUTO"))
