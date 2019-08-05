#' simulationFunctions.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' Dependencies:
library(adaptivetau)

#' Function list:
#'    stoich -- (not a function, but global definition)
#'    rateFunc()
#'    loadParameters()
#'    cellParameters()
#'    cellSSA()
#'    cellDivide()
#'    cellGenerate()
#'    simulateCells()
#'    simulateIC()


## 0. Basic stochastic simulation functions (stoichiometry, rate functions) ----
#' Define stoichiometry
#' Variable list:
#' x <- c("target.protein", "efflux.DNA", "efflux.mRNA", "efflux.protein", "drug", "target.dimer", "efflux.dimer", "alive")
#' nrow <- length(x)
#' Stoichiometry matrix
stoich <- matrix(
  c(
    #' TARGET dynamics
    +1, 0, 0, 0, 0, 0, 0, 0, # Constant target.protein translation
    #' EFFLUX dynamics
     0,+1, 0, 0, 0, 0, 0, 0, # efflux.DNA on
     0,-1, 0, 0, 0, 0, 0, 0, # efflux.DNA off
     0, 0,+1, 0, 0, 0, 0, 0, # efflux.mRNA transcription
     0, 0,-1, 0, 0, 0, 0, 0, # efflux.mRNA degradarion
     0, 0, 0,+1, 0, 0, 0, 0, # efflux.protein translation 
     0, 0, 0,-1, 0, 0, 0, 0, # efflux.protein degradation
    #' Drug dynamics
     0, 0, 0, 0,+1, 0, 0, 0, # Drug migration into cell
     0, 0, 0, 0,-1, 0, 0, 0, # Spontaneous loss of drug from cell
    #' target.dimer
    -1, 0, 0, 0,-1,+1, 0, 0, # Binding of target.protein and drug
    +1, 0, 0, 0,+1,-1, 0, 0, # Dissociation of target.dimer
    +1, 0, 0, 0, 0,-1, 0, 0, # Deacylation of target dimer
    #' efflux.dimer
     0, 0, 0,-1,-1, 0,+1, 0, # Binding of efflux.protein and drug
     0, 0, 0,+1,+1, 0,-1, 0, # Dissociation of target.dimer
     0, 0, 0,+1, 0, 0,-1, 0, # Removal/Catalysis of drug
     0, 0, 0,+1, 0, 0,-1, 0, # Deacylation of efflux dimer
    #' Cell death
     0, 0, 0, 0, 0, 0, 0,-1  # Cell death
  ), nrow = 8)

#' Define rate functions for our reaction scheme
#' 
#' \code{rateFunc}
#' 
#' @param x (named vector): a state vector of the form \code{x=c(target.protein=NN, efflux.dna=NN, efflux.mRNA=NN, efflux.protein=NN, drug=NN, target.dimer=NN, efflux.dimer=NN, alive={0,1})}
#' (NN = Nautral Number {0,1,2,3,...})
#' @param params (list): A list of parameters from \code{cell.parameters()}
#' @param t (numeric): the time, which is used to determine cell size
#' @return (vector): A numeric vector of transition rates
#' @examples rateFunc(x, params, t)
rateFunc <- function(x, params, t) {
  
  #' Define cell properties at time t
  cell.length <- with(params, cell.initialLength * 2 ^ (t / generationTime))
  cell.surfaceArea <- with(params, (pi * cell.diameter * cell.length) + cell.diameter / 2)
  cell.volume <- with(params, pi * (cell.diameter / 2) ^ 2 * cell.length)
  total.targets <- x[["target.protein"]] + x[["target.dimer"]]
  bound.fraction <- ifelse(total.targets == 0, 0, x[["target.dimer"]] / total.targets)
  
  #' Define parameters of death rate
  death.pars <- with(
    params,
    c(
      A = rate.maxGrowth * (rate.maxGrowth - rate.minGrowth) * (((1/MIC.fraction)^kappa) - 1) / (rate.maxGrowth * ((1/MIC.fraction)^kappa) - (rate.maxGrowth - rate.minGrowth)),
      B = rate.minGrowth * ((1/MIC.fraction)^kappa) / (rate.maxGrowth * ((1/MIC.fraction)^kappa) - (rate.maxGrowth - rate.minGrowth))
    )
  )
  #' Evaluate reaction rates
  rates <- with(
    params,
    c(
      #' TARGET reaction
      rate.target.protein.translation * x[["alive"]],
      #' EFFLUX reactions
      rate.efflux.DNA.switchOn * (1 - x[["efflux.DNA"]]) * x[["alive"]],
      rate.efflux.DNA.switchOff * x[["efflux.DNA"]] * x[["alive"]],
      rate.efflux.mRNA.transcription * x[["efflux.DNA"]] * x[["alive"]],
      rate.efflux.mRNA.deg * x[["efflux.mRNA"]] * x[["alive"]],
      rate.efflux.protein.translation * x[["efflux.mRNA"]] * x[["alive"]],
      rate.efflux.protein.deg * x[["efflux.protein"]] * x[["alive"]],
      #' Drug migration (Outside concentration is zero if t > drug.off time)
      ifelse(t < drug.off, drug.outsideConcentration, 0) * rate.drug.diffusion * cell.surfaceArea * x[["alive"]],
      (x[["drug"]] / cell.volume) * rate.drug.diffusion * cell.surfaceArea * x[["alive"]],
      #' TARGET dimer
      (rate.target.binding * x[["drug"]] * x[["target.protein"]] / cell.volume) * x[["alive"]],
      rate.dissociation * x[["target.dimer"]] * x[["alive"]],
      rate.deacylation * x[["target.dimer"]] * x[["alive"]],
      #' EFFLUX dimer
      (rate.efflux.binding * x[["drug"]] * x[["efflux.protein"]] / cell.volume) * x[["alive"]],
      rate.dissociation * x[["efflux.dimer"]] * x[["alive"]],
      rate.catalysis * x[["efflux.dimer"]] * x[["alive"]],
      rate.deacylation * x[["efflux.dimer"]] * x[["alive"]],
      #' Cell death (depends on fraction of bound targets, log(10) accounts for the fact that \psi are measured at the population level) [similar to Regoes 2004]
      log(10) * (death.pars[["A"]] * ((bound.fraction / MIC.fraction)^kappa)) / (((bound.fraction / MIC.fraction)^kappa) - death.pars[["B"]]) * x[["alive"]]
    )
  )
  return(rates)
}

## 1. Parameters: First loaded from file, and then modified for each individual cell ----

#' Load parameters from file.
#' Cell parameters are loaded from \code{"Params/cellParameters.R"}.
#' 
#' \code{load.parameters}
#' 
#' @param drug (string): Which drug parameters to use. This loads that parameter file.
#' This could be \code{"ciprofloxacin"}, \code{"rifampicin"}, or \code{"cefotaxime"}.
#' @return (list): Named list of parameters
#' @examples loadParameters(drug="rifampicin")
loadParameters <- function(drug = "ciprofloxacin") {
  #' Load cell parameters
  source("Params/cellParameters.R", local = TRUE)
  #' Load drug parameters
  source(paste("Params/", drug, ".R", sep = ""), local = TRUE)
  
  params <- list(
    #' TARGET parameters
    target.protein.total = drug.targets,
    rate.target.protein.translation = drug.targets / cell.generationTime, # Want average target density to be constant
    #' EFFLUX parameters
    rate.efflux.DNA.switchOn = efflux.DNA.switchOn,
    rate.efflux.DNA.switchOff = efflux.DNA.switchOff,
    rate.efflux.mRNA.transcription = efflux.mRNA.transcriptionRate,
    rate.efflux.mRNA.deg = efflux.mRNA.degRate,
    rate.efflux.protein.translation = efflux.protein.translationRate,
    rate.efflux.protein.deg = efflux.protein.degRate,
    #' Drug parameters
    drug.outsideConcentration = drug.outsideConcentration,
    rate.drug.diffusion = drug.diffusionConstant,
    MIC = drug.MIC,
    kappa = drug.kappa,
    rate.target.binding = drug.bindingRate,
    rate.efflux.binding = drug.bindingRate * efflux.binding.factor,
    rate.dissociation = dimer.dissociationRate,
    rate.deacylation = drug.deacylationRate,
    rate.catalysis = dimer.catalysisRate,
    #' Generation time, cell size, and growth rate params
    mean.generationTime = cell.generationTime,
    cell.diameter = cell.diameter,
    cell.initialLength = cell.initialLength,
    rate.maxGrowth = pop.maxGrowthRate,
    rate.minGrowth = pop.minGrowthRate,
    sd.generationTime = 0.05, # Percentage
    #' Mutant effects
    mutant.effect = mutant.effect,
    #' Mutant costs
    fitness.cost = fitness.cost,
    #' MIC fraction of bound targets
    MIC.fraction = MIC.fraction,
    #' Time to switch off the drug
    drug.off = Inf,
    #' Bias of efflux pumps to one pole during division [see Bergmiller et al. Science 356, 311 (2017)]
    efflux.binom.param = efflux.binom.param
  )
  return(params)
}

#' Parameters of individual cells, which are mutation-dependent or use random numbers.
#' 
#' \code{cell.parameters}
#' 
#' @param basic.params (list): A list of parameters from \code{loadParameters()}
#' @param mutation (string): Is our cell a mutant? What kind?
#' This can be \code{"WT"}, \code{"KO"}, \code{"REG-BURST"}, \code{"REG-ON"}, \code{"REG-OFF"}, \code{"STRUCT-BIND"}, or \code{"STRUCT-CAT"}.
#' @return (list): A named list in the same format of \code{loadParameters()}, but tailored to each cell.
#' @examples cellParameters(basic.params, mutation="REG-ON")
cellParameters <- function(basic.params, mutation = "WT") {
  #' Assign the basic model parameters
  cell.params <- basic.params
  
  #' Initialise the cost of mutant
  cost <- 0
  #' Now check the cell type, and apply conditions accordingly
  if (mutation == "KO") {
    #' Stop efflux gene activity
    cell.params$rate.efflux.DNA.switchOn <- 0.0
    #' Cost of knockout (is usually negative or zero, i.e. an advantage)
    cost <- -cell.params$fitness.cost
  }
  else if (mutation == "REG-ON") {
    #' Increase the switch-on rate of efflux DNA
    cell.params$rate.efflux.DNA.switchOn <- cell.params$rate.efflux.DNA.switchOn * cell.params$mutant.effect
    #' Cost of regulatory mutation
    cost <- cell.params$fitness.cost
  }
  else if (mutation == "REG-OFF") {
    #' Decrease the switch-off rate of efflux DNA
    cell.params$rate.efflux.DNA.switchOff <- cell.params$rate.efflux.DNA.switchOff / cell.params$mutant.effect
    #' Cost of regulatory mutation
    cost <- cell.params$fitness.cost
  }
  else if (mutation == "REG-BURST") {
    #' 1. Larger bursts of efflux mRNA
    # cell.params$rate.efflux.mRNA.transcription <- cell.params$rate.efflux.mRNA.transcription * cell.params$mutant.effect
    #' 2. Larger bursts of efflux Protein
    cell.params$rate.efflux.protein.translation <- cell.params$rate.efflux.protein.translation * cell.params$mutant.effect
    #' Cost of regulatory mutation
    cost <- cell.params$fitness.cost
  }
  else if (mutation == "STRUCT-BIND") {
    #' Increase binding rate of efflux to drug
    cell.params$rate.efflux.binding <- cell.params$rate.efflux.binding * cell.params$mutant.effect
    #' Cost of structural mutation
    cost <- cell.params$fitness.cost
  }
  else if (mutation == "STRUCT-CAT") {
    #' Increase catalysis rate of drug by efflux
    cell.params$rate.catalysis <- cell.params$rate.catalysis * cell.params$mutant.effect
    #' Cost of structural mutation
    cost <- cell.params$fitness.cost
  }
  else {#' mutation == "WT"
    cost <- 0
  }
  
  #' Fitness cost of mutations
  #' Increase average generation time for mutants, decrease for KO
  cell.params$mean.generationTime <- cell.params$mean.generationTime / (1.0 - cost)
  #' Modify all rates related to cell metabolism (here we have slowed down cell growth, and also translation)
  for (param in c("rate.target.protein.translation", "rate.efflux.protein.translation")) {
    cell.params[[param]] <- cell.params[[param]] * (1.0 - cost)
  }
  
  #' Generation time for each cell has a (uncorrelated) noise element (normally distributed)
  cell.params$generationTime <- rnorm(
    n = 1,
    mean = cell.params$mean.generationTime,
    sd = cell.params$sd.generationTime * cell.params$mean.generationTime
  )
  
  return(cell.params)
}


## 2. Simulation and division of a single cell ----

#' Execute the SSA for a single cell. Here we use the \code{adaptiveTau} package.
#' 
#' \code{cellSSA}
#' 
#' We can use two functions from the \code{adaptiveTau} package:
#' \code{ssa.exact()} uses the direct Gillespie simulation algorithm.
#' \code{ssa.adaptivetau()} uses the tau-leaping framework, and leads to faster (but apporximate, and maybe erroneous) simulations.
#' @param cell.params (list): A list of parameters for this cell from \code{cellParameters()}
#' @param initial.condition (named vector): A state vector of the form \code{x=c(target.protein=NN, efflux.dna=NN, efflux.mRNA=NN, efflux.protein=NN, drug=NN, target.dimer=NN, efflux.dimer=NN, alive={0,1})}
#' (NN = Nautral Number {0,1,2,3,...})
#' @param sim.t (numeric): How long to simulate for (in minutes). This is usually the generation time.
#' @return (dataframe): The timeseries of each molecule number as returned by the \code{adaptiveTau} algorithm.
#' @examples cellSSA(cell.params, initial.condition, sim.t=10)
cellSSA <- function(cell.params, initial.condition, sim.t) {
  #' Execute the SSA from t=0 until t=sim.t
  #cell.data <- ssa.exact(initial.condition, stoich, rateFunc, cell.params, sim.t)
  
  #' Execute the adaptiveTau simulation algorithm from t=0 to t=sim.t
  #' Here we use some error catching to reduce the chance of simulations failing
  cell.data <- tryCatch(
    #' First we try the algorithm with some default settings
    ssa.adaptivetau(initial.condition, stoich, rateFunc, cell.params, sim.t, tl.params = list(maxtau = 0.1, epsilon = 0.01)),
    #' If there is an error, we try to simulate with a greater accuracy.
    #' We do this recursively until we reach a threshold, at which poiint an error is thrown
    error = function(e) {
      epsilon <- 0.01
      #' Function to do the recursive calling
      robust <- function(epsilon) {
        epsilon <- epsilon / 2
        tryCatch(
          ssa.adaptivetau(initial.condition, stoich, rateFunc, cell.params, sim.t, tl.params = list(maxtau = 0.1, epsilon = epsilon)),
          error = function(e2) {
            if (epsilon < 1e-6) {
              stop(paste0("ERROR in adaptiveTau:\n", e2$message))
            }
            else {
              robust(epsilon)
            }
          }
        )
      }
      #' Execute this recursive function
      robust(epsilon)
    }
  )
  
  #' Convert to dataframe from matrix, and then return
  cell.data <- as.data.frame(cell.data)
  return(cell.data)
}

#' Split the within-cell molecules between the two daughter cells.
#' 
#' \code{cellDivide}
#' @param parent.cell.data (named vector): A state vector of the form \code{x=c(target.protein=NN, efflux.dna=NN, efflux.mRNA=NN, efflux.protein=NN, drug=NN, target.dimer=NN, efflux.dimer=NN, alive={0,1})}
#' (NN = Nautral Number {0,1,2,3,...})
#' @return (list): a list of two dataframes, each containing the protein numbers, ect. of each daughter in a single row
#' @examples cellDivide(parent.cell.data)
cellDivide <- function(parent.cell.data, parent.cell.params) {
  #' Cell division leads to binomially distributed attributes between the daughters (except DNA state, which is inherited from parent)
  #' The lines commented out below are options for evenly distributed proteins with and without 5% sd. 
  #' TARGET
  target.protein.1    <- rbinom(1, parent.cell.data[["target.protein"]], 0.5)
  #target.protein.1    <- round(rnorm(1, parent.cell.data[["target.protein"]] / 2, sd = 0.05 * parent.cell.data[["target.protein"]]))
  #target.protein.1    <- round(parent.cell.data[["target.protein"]] / 2)
  target.protein.2    <- parent.cell.data[["target.protein"]] - target.protein.1
  
  target.dimer.1      <- rbinom(1, parent.cell.data[["target.dimer"]], 0.5)
  #target.dimer.1    <- round(rnorm(1, parent.cell.data[["target.dimer"]] / 2, sd = 0.05 * parent.cell.data[["target.dimer"]]))
  #target.dimer.1    <- round(parent.cell.data[["target.dimer"]] / 2)
  target.dimer.2      <- parent.cell.data[["target.dimer"]] - target.dimer.1
  
  #' EFFLUX
  efflux.DNA.1        <- parent.cell.data[["efflux.DNA"]] # DNA states are directly inherited from parent
  efflux.DNA.2        <- parent.cell.data[["efflux.DNA"]]
  
  efflux.mRNA.1       <- rbinom(1, parent.cell.data[["efflux.mRNA"]], 0.5)
  efflux.mRNA.2       <- parent.cell.data[["efflux.mRNA"]] - efflux.mRNA.1
  
  
  efflux.protein.1    <- rbinom(1, parent.cell.data[["efflux.protein"]], parent.cell.params$efflux.binom.param)
  #efflux.protein.1    <- rbinom(1, parent.cell.data[["efflux.protein"]], 0.5) # Originally we assumed the binomial parameter is 0.5, but we now add this as an extra parameter
  #efflux.protein.1    <- round(rnorm(1, parent.cell.data[["efflux.protein"]] / 2, sd = 0.05 * parent.cell.data[["efflux.protein"]]))
  #efflux.protein.1    <- round(parent.cell.data[["efflux.protein"]] / 2)
  efflux.protein.2    <- parent.cell.data[["efflux.protein"]] - efflux.protein.1
  
  efflux.dimer.1      <- rbinom(1, parent.cell.data[["efflux.dimer"]], parent.cell.params$efflux.binom.param)
  #efflux.dimer.1      <- rbinom(1, parent.cell.data[["efflux.dimer"]], 0.5) # Originally we assumed the binomial parameter is 0.5, but we now add this as an extra parameter
  #efflux.dimer.1      <- round(rnorm(1, parent.cell.data[["efflux.dimer"]] / 2, sd = 0.05 * parent.cell.data[["efflux.dimer"]]))
  #efflux.dimer.1      <- round(parent.cell.data[["efflux.dimer"]] / 2)
  efflux.dimer.2      <- parent.cell.data[["efflux.dimer"]] - efflux.dimer.1
  
  #' Drug
  drug.1              <- rbinom(1, parent.cell.data[["drug"]], 0.5)
  drug.2              <- parent.cell.data[["drug"]] - drug.1
  
  #' Create dataframes for each daughters' state variables, and return as a list
  daughter.1.data <- data.frame(
    time = parent.cell.data[["time"]],
    target.protein = target.protein.1,
    efflux.DNA = efflux.DNA.1,
    efflux.mRNA = efflux.mRNA.1,
    efflux.protein = efflux.protein.1,
    drug = drug.1,
    target.dimer = target.dimer.1,
    efflux.dimer = efflux.dimer.1,
    alive = 1
  )
  daughter.2.data <- data.frame(
    time = parent.cell.data[["time"]],
    target.protein = target.protein.2,
    efflux.DNA = efflux.DNA.2,
    efflux.mRNA = efflux.mRNA.2,
    efflux.protein = efflux.protein.2,
    drug = drug.2,
    target.dimer = target.dimer.2,
    efflux.dimer = efflux.dimer.2,
    alive = 1
  )
  daughter.cells <- list(daughter.1.data, daughter.2.data)
  return(daughter.cells)
}


## 3. Generate a blank cell, useful for initial initial conditions ----

#' Initiate a blank cell with no molecules.
#' Useful for generating initial conditions.
#' 
#' \code{cellGenerate()}
#' @params (list): A list of parameters for this cell from \code{cellParameters()}
#' @return (dataframe): A state vector of the form \code{x=c(target.protein=X, efflux.dna=0, efflux.mRNA=0, efflux.protein=0, drug=NN, target.dimer=0, efflux.dimer=0, alive=1)}
#' (X = number of targets from cell.params)
#' @examples cellGenerate(cell.params)
cellGenerate <- function(cell.params) {
  
  cell.data <- data.frame(
    target.protein = cell.params$target.protein.total,
    efflux.DNA = 0,
    efflux.mRNA = 0,
    efflux.protein = 0,
    drug = 0,
    target.dimer = 0,
    efflux.dimer = 0,
    alive = 1
  )
  return(cell.data)
}

## 4. Master simulation code ----

#' Run simulations in which cells are simulated chronologically, and then divided.
#' This runs until the time exceeds \code{t.max}.
#' 
#' \code{simulateCells()}
#' @param basic.params (list): A list of parameters from \code{loadParameters()}
#' @param initial.cell.data (dataframe): Contains info about protein numbers, etc. in the initial cells (not time!)
#' @param initial.cell.types (string vector): A list of whether cells are WT or mutant (Reg/Struct)
#' @param t.max (numeric): Defines how long the simulation should run for (in minutes).
#' @param max.num.cells (numeric): How many cells are allowed at one time before we stop the simulations?
#' @param single.cell (Bool): If TRUE, we only track a single cell, rather than the dynamic population.
#' @param store.info (Bool): If TRUE, we store the within-cell data.
#' @return (list): A list of the following objects:
#' @return cell.params.list (list): List of list of cell parameters
#' @return cell.data.list (list): List of dataframes of each cell's within-cell variables (protein number, drug dynamics, alive/dead status, etc. over time)
#' @return family.tree (dataframe): A dataframe of each cell's fundamental properties (ID, parent, birth time, end of life, did it die, is it a mutant...)
#' @examples simulateCells(basic.params, initial.cell.data, initial.cell.types=c("WT","WT","REG-ON"), t.max=100)
simulateCells <- function(basic.params,
                          initial.cell.data,
                          initial.cell.types,
                          t.max,
                          max.num.cells = 1000,
                          single.cell = FALSE,
                          store.info = TRUE
) {
  
  #' Initialisation
  #' Start time
  t.0 <- 0.0
  
  #' Identifiers of the initial cells
  cell.id <- as.character(seq_along(initial.cell.types))
  
  # Cell parameter list
  #' Each cell has its own parameters, which we store as a list
  #' The name of the list element is the cell.id
  #' We use cell.parameters() to generate each parameter set for our initial cells
  cell.params.list <- lapply(initial.cell.types, function(mutation) cellParameters(basic.params, mutation))
  names(cell.params.list) <- cell.id
  
  # Cell data list
  #' A list of dataframes, each describing the dynamics within each cell over its lifetime.
  #' First we need to add a time column to our initial.cell.data dataframe
  initial.cell.data$time <- t.0
  #' The name of the list element is the cell.id
  #' We split the initial.cell.data dataframe into single rows, and create a list
  cell.data.list <- lapply(cell.id, function(i) data.frame(initial.cell.data[i,], row.names = NULL))
  names(cell.data.list) <- cell.id
  
  # Family tree
  #' A dataframe containing cell properties
  family.tree <- data.frame(
    cell.id = cell.id,
    parent.id = "0",
    first.time = t.0,
    last.time = NA,
    alive = NA,
    mutation = unlist(initial.cell.types),
    generation.time = unlist(
      lapply(cell.params.list, function(cell)
        cell$generationTime),
      use.names = FALSE
    ),
    simulated = 0,
    
    stringsAsFactors = FALSE
  )
  
  #' We have a cutoff which stops our simulations exploding.
  #' If more than max.num.cells cells are alive, we stop
  num.cells <- nrow(family.tree)
  
  #' Loop over time and let our cells grow! (or die)
  min.t <- t.0
  while (min.t <= t.max) {
    #' Choose cell to simulate
    #' First row in which we have minimum time and has not yet been simulated
    this.cell.id <- as.character(subset(family.tree, subset = (first.time == min.t & simulated == 0))[1, "cell.id"])
    
    #' Simulate cell
    #' Get parameters
    this.cell.params <- cell.params.list[[this.cell.id]]
    #' Get initial condition
    this.cell.data.initial <- cell.data.list[[this.cell.id]]
    #' Turn into vector (without time column)
    this.cell.data.initial <- unlist(this.cell.data.initial[nrow(this.cell.data.initial), names(this.cell.data.initial) != "time"])
    #' Determine time to simulate for:
    #' First determine generation time, then time to end of simulation.
    #' We simulate for the minimum of the two.
    this.cell.generation.time <- this.cell.params$generationTime
    time.to.end <- t.max - min.t
    sim.t <- min(this.cell.generation.time, time.to.end)
    #' Does the drug switch off in this cells lifetime?
    this.cell.params$drug.off <- basic.params$drug.off - min.t
    #' Run simulation
    this.cell.data <- cellSSA(this.cell.params, this.cell.data.initial, sim.t)
    #' If cell died, remove last row
    if (this.cell.data[nrow(this.cell.data), "alive"] == 0) {
      this.cell.data <- this.cell.data[c(1:(nrow(this.cell.data) - 1)), ]
    }
    #' Update time in this dataframe
    this.cell.data$time <- this.cell.data$time + min.t
    #' Store this simulation data in cell.data.list (if not flagged out)
    if (store.info) cell.data.list[[this.cell.id]] <- this.cell.data
    
    #' Update the information about this cell to the family tree
    #' First get the information
    this.cell.information <- family.tree[family.tree$cell.id == this.cell.id, ]
    #' Update final time seen by this cell
    this.cell.information$last.time <- this.cell.data[nrow(this.cell.data), "time"]
    #' Update whether this cell died or not
    this.cell.information$alive <- this.cell.data[nrow(this.cell.data), "alive"]
    #' Update that this cell has been simulated
    this.cell.information$simulated <- 1
    #' Finally, reassign information to family.tree dataframe
    family.tree[family.tree$cell.id == this.cell.id, ] <- this.cell.information
    
    
    #' Divide cells
    #' First test if we should divide:
    #' We should only divide if the cell is alive,
    #' and the entire generation time was simulated.
    #' We only keep one daughter if [single.cell==TRUE]
    if (this.cell.information$alive == 1 & sim.t == this.cell.generation.time) {
      #' Extract the parent cell data
      parent.cell.data <- unlist(this.cell.data[nrow(this.cell.data),])
      parent.cell.information <- this.cell.information
      #' Perform division and get data for two daughters
      daughters.data.list <- cellDivide(parent.cell.data, this.cell.params)
      #' Drop a daughter if single.cell==TRUE
      if (single.cell) daughters.data.list <- daughters.data.list[1]
      #' Define ID values for these daughters
      max.id <- max(as.numeric(family.tree$cell.id))
      daughters.id <- as.character(max.id + c(1:length(daughters.data.list)))
      names(daughters.data.list) <- daughters.id
      #' Add these daughters to cell.data.list
      cell.data.list <- c(cell.data.list, daughters.data.list)
      #' Draw parameters for these cells, name, and add to cell.params.list
      daughters.params.list <- lapply(daughters.data.list, function(i) cellParameters(basic.params, mutation = parent.cell.information$mutation))
      names(daughters.params.list) <- daughters.id
      cell.params.list <- c(cell.params.list, daughters.params.list)
      #' Add family tree entries for the daughter cells
      daughter.information <- data.frame(
        cell.id = daughters.id,
        parent.id = this.cell.id,
        first.time = parent.cell.information$last.time,
        last.time = NA,
        alive = NA,
        mutation = parent.cell.information$mutation,
        generation.time = unlist(
          lapply(daughters.params.list, function(cell)
            cell$generationTime),
          use.names = FALSE
        ),
        simulated = 0,
        
        stringsAsFactors = FALSE
      )
      family.tree <- rbind(family.tree, daughter.information)
    }
    
    #' Simulation time is based on the earliest birth time of cells not yet simulated
    #' Returns Inf if empty to break the loop
    min.t <- suppressWarnings(min(family.tree[family.tree$simulated == 0, "first.time"]))
    #' Total number of cells alive at this time
    num.cells <- sum(family.tree$simulated == 0)
    if (num.cells >= max.num.cells) t.max <- max(family.tree[, "first.time"])
  }
  
  #' Output results
  out <- list(
    family.tree = family.tree,
    cell.data = cell.data.list,
    cell.params = cell.params.list
  )
  return(out)
}


#' Simulate or load the initial conditions in the absence of drug
#' 
#' \code{simulateCells()}
#' @param basic.params (list): A list of parameters from \code{loadParameters()}
#' @param num.cells (named numeric): How many cells do we want to initialise? The names of the vector are the cell types.
#' @param dir (char): The filename from where we try to load initial conditions, if they are already stored.
#' @param t.max (numeric): Defines how long the simulation should run for (in minutes).
#' @return (list): A list of the following objects:
#' @return cell.params.list (list): List of list of cell parameters
#' @return cell.data.list (list): List of dataframes of each cell's within-cell variables (protein number, drug dynamics, alive/dead status, etc. over time)
#' @return family.tree (dataframe): A dataframe of each cell's fundamental properties (ID, parent, birth time, end of life, did it die, is it a mutant...)
#' @examples simulateCells(basic.params, initial.cell.data, initial.cell.types=c("WT","WT","REG-ON"), t.max=100)
simulateIC <- function(basic.params, num.cells, dir="false", t.max=200) {
  #' Create vector of initial cell types
  initial.cell.types <- rep(names(num.cells), times = num.cells)
  
  #' Identify file names (if possible)
  filename <- paste("initialCondition_", names(num.cells), ".dat", sep = "")
  files <- paste(dir, filename, sep = "/")
  names(files) <- names(num.cells)
  
  #' Dataframe for all ICs
  initial.cell.data <- cellGenerate(basic.params)[0,]
  #' If file exists, then load the file. Else generate cells
  for (cell.type in names(num.cells)) {
    if (num.cells[[cell.type]] > 0) {
      #' Does the initial condition file exist
      if (file.exists(files[[cell.type]])) {
        #' Load file
        IC.data <- read.csv(files[[cell.type]])
        #' Sample the ICs
        IC.data <- IC.data[sample(nrow(IC.data), num.cells[[cell.type]]), ] 
      }
      else{
        #' Generate initial cells using single-cell simulations.
        #' Create parameter set without drug
        basic.params$drug.outsideConcentration <- 0
        #' Cell to initialise
        initial.cell <- cellGenerate(basic.params)
        #' Generate cells
        IC.data <- DFapply(seq_len(num.cells[[cell.type]]), function(i) {
          sim.result <- simulateCells(basic.params, initial.cell, cell.type, t.max = t.max, single.cell = TRUE)
          #' Return only last line
          return(sim.result$cell.data[[length(sim.result$cell.data)]][1,])
        })
        #' Drop time column
        IC.data <- IC.data[, names(IC.data) != "time"]
      }
      #' Store
      initial.cell.data <- rbind(initial.cell.data, IC.data, make.row.names = FALSE)
    }
  }
  
  return(initial.cell.data)
}
