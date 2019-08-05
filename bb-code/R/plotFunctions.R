#' plotFunctions.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' Dependencies:
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(cowplot)

#' Function list:
#' - palette_OkabeIto
#' - gradientPalette()
#' - label_panels()
#' - sciFormat()
#' - logTicks
#' - cellColours()
#' - plotSave()
#' - plotPopSize()
#' - plotPopSizeMulti()
#' - plotWithinCell()

#' Nice colour scheme defintion
palette_OkabeIto <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#' Gradient colour scheme
gradientPalette <- function(n, palette = "Blues") {
    colours <- brewer.pal(9, palette)[5:9]
    colours <- colorRampPalette(colours)(n)
    return(colours)
}

#' Function to add A, B, C, etc. to a facet label
#' @examples facet_wrap(~facet, labeller = label_panels)
make_labelstring <- function(mypanels) {
    mylabels <- sapply(mypanels, function(x) paste(LETTERS[which(mypanels == x)], x, sep = ": "))
    return(mylabels)
}
label_panels <- as_labeller(make_labelstring)

#' Function for scientific notation on y-axis
#'
#' \code{sciFormat()}
#' @param l (string): Number to be converted to expoenential format
#' @return (string)
#' @examples sciFormat("1e-1")
sciFormat <- function(l) {
    l <- format(l, scientific = TRUE)
    ## Just plot 10^x
    l <- gsub("e","10^", gsub("^(.*)e", "e", l))
    # Return this as an expression
    parse(text = l)
}

#' Function to return major and minor logarithmic tick points
#' 
#' \code{logTicks()}
#' @param minor (bool): Do we return minor ticks as well as major?
#' @return (list): A list of two vectors (major and minor tick values)
#' @examples logTicks(minor = TRUE)
logTicks <- function(minor = TRUE) {
    breaks <- 10^c(-4:9)
    if (minor) {
        minor.breaks <- unlist(lapply(breaks, function(lower) lower * c(2:9) ))
    }
    else minor.breaks <- NULL
    return(list(major = breaks, minor = minor.breaks))
}


#' Assign colours to each cell type
cellColours <- function() {
    #c(KO = "darkmagenta", WT = "lightblue", `REG-ON` = "orange", `REG-OFF` = "orange", `REG-BURST` = "orange", `STRUCT-BIND` = "green", `STRUCT-CAT` = "green")
    colours <- palette_OkabeIto[1:7]
    #names(colours) <- c("WT", "STRUCT-CAT", "STRUCT-BIND", "REG-ON", "REG-OFF", "REG-BURST", "KO")
    names(colours) <- c("KO", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT", "REG-ON", "REG-OFF", "WT")
    return(colours)
}

#' My saving function to store figures with nice dimensions
#' 
#' \code{plotSave}
#' 
#' @param filename (string): name to be assigned to the stored figure.
#' @param plot (gg-object): the figure to be saved.
#' @param width (numeric): The width of the fiugure in millimeters.
#' @param height.mult (numeric): Ratio of height to width. Defaults to golden ratio.
#' @param ... : Options for ggsave()
#' @return NULL
#' @examples plotSave("./foo.pdf", plot, width=20, height.mult=2)
plotSave <- function(filename = "~/Desktop/foo.pdf", plot, width, height.mult = 2/(1 + sqrt(5)), ...){
    ggsave(filename = filename, plot = plot, width = width, height = width*height.mult, units = "mm", ...)
}

#' Plot the population size curves over time
#' \code{plotPopSize}
#' 
#' @param pop.size (dataframe): Size of population over time, broken up into cell groups
#' @param stack.area (boolean): If true, we stack each population on top of each other,
#' otherwise we just plot absolute number as well as the total
#' @return A graphical object showing a skyline plot of cell numbers
#' @examples plotPopSize(pop.size)
plotPopSize <- function(pop.size, stack.area=FALSE) {
    # Which variables are we plotting
    not.plot <- c("time", "sim.label")
    if (stack.area) not.plot <- c(not.plot, "Total")
    measure.vars <- names(pop.size)[!names(pop.size) %in% not.plot]
    
    # Assign colours for each cell type
    colour.def <- c(cellColours(), Total = "black")
    colours <- colour.def[measure.vars]
    
    # Melt dataframe for easy plotting
    pop.size.melt <- melt(
        data = pop.size,
        id.vars = "time",
        measure.vars = measure.vars,
        variable.name = "genotype",
        value.name = "number"
    )
    # Generate the skyline plot
    plot <- ggplot(pop.size.melt, aes(x = time, y = number, colour = genotype, fill = genotype)) +
        geom_step(position = ifelse(stack.area, "stack", "identity")) +
        theme_bw() +
        scale_color_manual(values = colours) +
        scale_fill_manual(values = colours)
    return(plot)
}

#' Plot multiple population size curves over time
#' 
#' \code{plotPopSizeMulti}
#' 
#' @param pop.size.list (list of dataframes): Size of population over time, broken up into cell groups
#' @param stack.area (boolean): If true, we stack each population on top of each other,
#' otherwise we just plot absolute number as well as the total
#' @return A graphical object showing a skyline plot of cell numbers
#' @examples plotPopSizeMulti(pop.size.list)
plotPopSizeMulti <- function(pop.size.list, stack.area=FALSE) {
    # Add simulation label to each dataframe in the list
    for (i in seq_along(pop.size.list)) pop.size.list[[i]]$sim.label <- i
    # Combine
    pop.size <- do.call(rbind, pop.size.list)
    
    # Which variables are we plotting
    not.plot <- "time"
    if (stack.area) not.plot <- c(not.plot, "Total")
    measure.vars <- names(pop.size)[!names(pop.size) %in% not.plot]
    
    # Assign colours for each cell type
    colour.def <- c(cellColours(), Total = "black")
    colours <- colour.def[measure.vars]
    
    # Melt dataframe for easy plotting
    pop.size.melt <- melt(
        data = pop.size,
        id.vars = c("time", "sim.label"),
        measure.vars = measure.vars,
        variable.name = "genotype",
        value.name = "number"
    )
    # Generate the skyline plot
    plot <- ggplot(pop.size.melt, aes(x = time, y = number, colour = genotype)) +
        geom_step(position = ifelse(stack.area, "stack", "identity") ) +
        facet_wrap(~sim.label) +
        theme_bw() +
        scale_color_manual(values = colours)
    return(plot)
}

#' Plot the timeseries of within-cell components
#' @param ssa.data (dataframe): Data from the SSA, with added generation labels, concentrations, and fraction of bound targets
#' @param ode.data (dataframe): Data from the ODE model, with added generation labels, concentrations, and fraction of bound targets
#' @param measure.vars (vector): Which quantities do we want to plot. If named, then we rename the factors
#' @value ggplot
plotWithinCell <- function(ssa.data, ode.data, measure.vars, ncol=NULL) {
    #' Melt the SSA and ODE data, ready for plotting
    meltData <- function(data, measure.vars) melt(data = data, id.vars = c("time", "generation"), measure.vars = measure.vars, variable.name = "molecule", value.name = "number")
    ssa.data <- meltData(ssa.data, measure.vars)
    ode.data <- meltData(ode.data, measure.vars)
    
    #' If named, then reset the factors
    if (!is.null(names(measure.vars))) {
        levels(ssa.data$molecule) <- names(measure.vars)
        levels(ode.data$molecule) <- names(measure.vars)
    }
    
    #' Plot
    plot <- ggplot(ssa.data, aes(x = time, y = number, colour = generation)) +
        geom_path(data = ode.data, aes(x = time, y = number), colour = palette_OkabeIto[8], inherit.aes = FALSE) +
        geom_step() +
        #ifelse(is.null(ncol), facet_wrap(~molecule, scales = "free_y"), facet_wrap(~molecule, scales = "free_y", ncol = ncol)) +
        scale_color_manual(values = rep_len(c(palette_OkabeIto[5], palette_OkabeIto[6]), length(levels(ssa.data$generation)))) +
        theme_bw() +
        theme(legend.position = "None",
              strip.background = element_blank(),
              strip.text = element_text(hjust = 0),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        labs(x = "time (minutes)", y = NULL)
    
    if (is.null(ncol)) plot <- plot + facet_wrap(~molecule, scales = "free_y")
    else plot <- plot + facet_wrap(~molecule, scales = "free_y", ncol = ncol)
        
    return(plot)
}
