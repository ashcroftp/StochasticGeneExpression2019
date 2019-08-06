#' get-params.R
#' Author: Peter Ashcroft, ETH Zurich

#' A script to calculate MIC fraction of multiple cell types based on the mean-field equations

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/meanFieldModel.R")
source("R/dataAnalysisFunctions.R")
source("R/plotFunctions.R")
library(GGally)

# The cell parameters we consider ----
new.cell.params <- expand.grid(
    list(
        binding = c(1e4, 1e5, 1e6),
        catalysis = c(1e-4, 1e-2, 1e1),
        diffusion = c(1e-11, 1e-10, 1e-9),
        targets = c(10, 100, 1000),
        efflux.pumps = c(10, 100, 1000),
        drug.molecules = c(10, 100, 1000)
    ),
    KEEP.OUT.ATTRS = F)

# Some constants used to generate the data
#' Define some initial cell types and basic parameters
cell.type <- "WT"
initial.cell <- c(1)
names(initial.cell) <- cell.type
basic.params <- loadParameters("ciprofloxacin")
t.max <- 400

# Simulate ----
#' Modify for each new parameter set, then integrate the ODE and output values
params.output <- parallel.DFapply(seq_len(nrow(new.cell.params)), function(i) {
    
    #' Modify the parameters based on the new cell properties
    this.params <- basic.params
    #' Binding rates
    this.params$rate.target.binding <- new.cell.params[i, "binding"] * 60.0*1.0e15/6.022e23
    this.params$rate.efflux.binding <- this.params$rate.target.binding
    #' Catalysis rate
    this.params$rate.catalysis <- new.cell.params[i, "catalysis"] * 60.0
    #' Diffusion rate
    this.params$rate.drug.diffusion <- new.cell.params[i, "diffusion"] * 60.0*1.0e6
    #' Number of targets
    this.params$target.protein.total <- round(new.cell.params[i, "targets"] * (2/3))
    this.params$rate.target.protein.translation <- this.params$target.protein.total / this.params$mean.generationTime
    #' Number of efflux pumps
    this.params$rate.efflux.mRNA.transcription <- new.cell.params[i, "efflux.pumps"] / 15.616
    
    #' The drug concentration is determined based on numerical methods.
    #' We integrate the ODEs with different external concentrations until the cell contains
    #' an average of new.cell.params$drug.molecules[i] molecules of drug.
    #' We use a secant method to find the external concentration

    #' Simulate the ODEs
    set.seed(11)
    initial.cell.data <- simulateIC(this.params, initial.cell)
    #' Function to integrate the ODEs and return the number of drug molecules per cell
    getDrug <- function(drug) {
        this.params$drug.outsideConcentration <- drug
        ode.sol <- integrateDiluteEqns(t.max, params = this.params, IC = initial.cell.data)
        sol <- ode.sol[ode.sol$time > t.max - 2 * this.params$mean.generationTime, ]
        drug.inside <- mean(sol$drug + sol$target.dimer + sol$efflux.dimer)
        return(drug.inside)
    }    
    #' Initial drug concentration guesses
    x <- c(0.01, 1000 * new.cell.params[i, "drug.molecules"])
    y <- (c(getDrug(x[1]), getDrug(x[2])) - new.cell.params[i, "drug.molecules"]) / new.cell.params[i, "drug.molecules"]
    epsilon <- 1
    j <- 2
    while (epsilon > 0.001) {
        j <- j + 1
        x[j] <- x[j - 1] - y[j - 1] * ((x[j - 1] - x[j - 2])/(y[j - 1] - y[j - 2]))
        y[j] <- (getDrug(x[j]) - new.cell.params[i, "drug.molecules"]) / new.cell.params[i, "drug.molecules"]
        epsilon <- abs(y[j])
    }
    #' We assign the solved drug concentration as the MIC for this particular cell
    this.params$MIC <- x[length(x)]
    rm(j,x,y)
    
    #' We then integrate the ODEs at this MIC value and then extract the equilibrium values, from which we compute the MIC fraction
    this.params$drug.outsideConcentration <- this.params$MIC
    ode.sol <- integrateDiluteEqns(t.max, params = this.params, IC = initial.cell.data)
    #' For the equilibrium values, we average over the last two generations only
    sol <- ode.sol[ode.sol$time > t.max - 2 * this.params$mean.generationTime, ]
    #' Extract MIC fraction
    this.params$MIC.fraction <- mean(sol$target.dimer / (sol$target.protein + sol$target.dimer))
    
    #' Now create a dataframe with the cell's prescribed values and the modified parameters
    out <- cbind(
        new.cell.params[i, ],
        data.frame(
            #' Binding rates
            rate.target.binding = this.params$rate.target.binding,
            rate.efflux.binding = this.params$rate.efflux.binding,
            #' Catalysis rate
            rate.catalysis = this.params$rate.catalysis,
            #' Diffusion rate
            rate.drug.diffusion = this.params$rate.drug.diffusion,
            #' Number of targets
            target.protein.total = this.params$target.protein.total,
            rate.target.protein.translation = this.params$rate.target.protein.translation,
            #' Number of efflux pumps
            rate.efflux.mRNA.transcription = this.params$rate.efflux.mRNA.transcription,
            #' Drug MIC concentration
            MIC = this.params$MIC,
            #' MIC fraction
            MIC.fraction = this.params$MIC.fraction
        )
    )
    return(out)
})


write.csv(params.output, "./Systematic/zz-data/params.output", row.names = F)


# Perform linear regression ----
#' First inspect the results vs variables
melted <- melt(data = cbind(params.output, data.frame(colour = params.output$MIC.fraction)), id.vars = c("MIC", "MIC.fraction","colour"), measure.vars = c("binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules"))
ggplot(melted, aes(x = value, y = MIC, colour = colour)) +
    geom_smooth(method = "lm", se = F, colour = palette_OkabeIto[8]) +
    geom_point() +
    scale_color_viridis_c(name = "MIC fraction", limits = c(0,1)) +#(trans = "log10") +
    facet_grid(~variable, scales = "free_x", labeller = label_panels) +
    scale_x_log10(labels = sciFormat) +
    scale_y_log10(labels = sciFormat) + 
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

#' Try pairwise comparisons
dat <- params.output[, c("MIC", "MIC.fraction", "binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules")]
for (j in c("binding", "catalysis", "diffusion", "targets", "efflux.pumps", "drug.molecules", "MIC")) {
    dat[, j] <- log10(dat[, j])
}
rm(j)
ggpairs(dat) +
    theme_bw()  +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )


mod <- lm(MIC ~ binding + catalysis + diffusion + targets + efflux.pumps + drug.molecules, data = params.output)
pred <- predict.lm(mod)
dat <- data.frame(MIC = params.output$MIC, pred = pred)

ggplot(dat, aes(x = MIC, y = pred)) + geom_point() + geom_abline(intercept = 0, slope = 1)

# # Plot ----
# variables <- c(
#     `number of free drug molecules` = "drug",
#     `number of unbound target proteins` = "target.protein",
#     `number of bound target proteins` = "target.dimer",
#     `fraction of bound target proteins` = "fraction.bound.targets",
#     `efflux DNA state` = "efflux.DNA",
#     `number of efflux mRNA` = "efflux.mRNA",
#     `number of efflux proteins` = "efflux.protein",
#     `number of bound efflux proteins` = "efflux.dimer"
# )
# melted <- DFapply(names(ode.sol), function(type) {
#     df <- melt(data = ode.sol[[type]], id.vars = c("time", "generation"), measure.vars = variables, variable.name = "molecule", value.name = "number")
#     df$type <- factor(type, levels = names(ode.sol))
#     levels(df$molecule) <- names(variables)
#     return(df)
# })
# ggplot(melted, aes(x = time, y = number, colour = type)) +
#     geom_step() +
#     facet_wrap(~molecule, scales = "free_y", nrow = 4, labeller = label_panels) +
#     theme_bw() +
#     theme(legend.position = "None",
#           strip.background = element_blank(),
#           strip.text = element_text(hjust = 0),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()) +
#     labs(x = "time (minutes)", y = NULL)
