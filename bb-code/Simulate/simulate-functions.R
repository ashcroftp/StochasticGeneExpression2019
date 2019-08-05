#' simulate-functions.R
#' Author: Peter Ashcroft, ETH Zurich

#' Functions used when submitting jobs to EULER etc.

#' Dependencies:
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/checks.R")
source("Simulate/filenames.R")
library(data.table)

#' Function list:
#' - getIndex()
#' - getIndices()
#' - getExpectedFiles()
#' - getMissingFiles()
#' - getMissingSims()
#' - modifyParams()
#' - runSim()
#' - runSurvival()
#' - runGrowth()


#' Compute the combined parameter--node index (inverse of getIndices())
getIndex <- function(param.index, node.index, num.params, num.nodes) {
    index <- (param.index - 1) * num.nodes + node.index
    return(index)
}


#' Separate the parameter and node indices from the combined index (inverse of getIndex())
getIndices <- function(index, num.params, num.nodes) {
    param.index <- ceiling(index / num.nodes)
    node.index <- index %% num.nodes
    if (node.index == 0) node.index <- num.nodes
    return(c(param = param.index, node = node.index))
}


#' Which output files do we expect to find in the data directory
getExpectedFiles <- function(data.directory, file.id, num.params, num.nodes, file.type = "output") {
    df <- expand.grid(node.index = seq_len(num.nodes), param.index = seq_len(num.params))
    df <- df[, c("param.index", "node.index")]
    #' Add the EULER index
    df$index <- getIndex(df$param.index, df$node.index, num.params, num.nodes)
    #' Add the filename
    df$file <- sapply(seq_len(nrow(df)), function(i) {
        if (file.type == "output") {
            output.file(data.directory, file.id, df[i, "param.index"], df[i, "node.index"])
        } else if (file.type == "seed") {
            seed.file(data.directory, file.id, df[i, "param.index"], df[i, "node.index"])
        }
    }, USE.NAMES = F)
    return(df)
}


#' Which output files are missing from the data directory
getMissingFiles <- function(data.directory, file.id, num.params, num.nodes) {
    #' List the expected output files
    expected.files <- getExpectedFiles(data.directory, file.id, num.params, num.nodes)
    
    #' List the stored output files
    output.files <- paste0(data.directory, list.files(data.directory, pattern = "*_output_*"))
    
    #' Which files are missing
    expected.files$missing <- sapply(seq_len(nrow(expected.files)), function(i) !expected.files[i, "file"] %in% output.files)
    return(expected.files[expected.files$missing, c("param.index", "node.index", "index")])
}


#' Which simulations are missinbg from the output files?
getMissingSims <- function(data.directory, file.id, num.params, num.nodes, sims.per.node) {
    #' We recover the vector of completed simulation labels from the seed file
    #' We then use this to determine which simulations are missing from each file, and re-run these accordingly
    
    #' List the expected seed files
    seed.files <- getExpectedFiles(data.directory, file.id, num.params, num.nodes, file.type = "seed")
    
    #' Read each seed file in turn, and if simulations are missing we return a list with elements param.index, node.index, and a vector = missing sim.labels
    missing.labels <- sapply(seq_len(nrow(seed.files)), function(i) {
        #' Compute the expected sim.labels
        expected.labels <- (seed.files[i, "node.index"] - 1) * sims.per.node + seq_len(sims.per.node)
        #' Read the file and extract the sim.label column
        sim.labels <- fread(seed.files[i, "file"], sep = ",", select = c("sim.label"))$sim.label
        #' Which labels are missing
        missing.labels <- expected.labels[!expected.labels %in% sim.labels]
        if (length(missing.labels) > 0) {
            return(list(indices = c(param.index = seed.files[i, "param.index"], node.index = seed.files[i, "node.index"]), sim.labels = missing.labels))
        } else {
            return(NULL)
        }
    }, simplify = F)
    #' Add a name to the vector
    names(missing.labels) <- seed.files$index
    #' Remove NULL elements
    missing.labels[sapply(missing.labels, is.null)] <- NULL
    return(missing.labels)
}


#' Incorporate parameter changes into the local parameters
modifyParams <- function(basic.params, var.params, var.param.index){
    for (param in names(var.params)) {
        #' Is this parameter a multiplicative factor?
        if (param != strsplit(param, ".MULT")[[1]]) {
            if (param == "drug.outsideConcentration.MULT") {
                #' If drug concentration, then we work in units of MIC
                basic.params[strsplit(param, ".MULT")[[1]]] <- basic.params$MIC * var.params[var.param.index, param]
            } else {
                #' Otherwise, we just multiply the original value
                basic.params[strsplit(param, ".MULT")[[1]]] <- basic.params[[strsplit(param, ".MULT")[[1]]]] * var.params[var.param.index, param]
            }
        } else {
            #' If not multiplicative, then just reset the parameter value
            basic.params[param] <- var.params[var.param.index, param]
        }
    }
    rm(param)
    return(basic.params)
}


#' Repeatedly call the simulation algorithm
runSim <- function(data.directory, file.id, cell.type, output, num.params, num.nodes, sims.per.node, index, local.parallel = FALSE, sim.labels = NULL) {
    #' Convert the unique index into parameter and node indices
    indices <- getIndices(index, num.params, num.nodes)
    
    #' Load parameters from local file, and modify parameters based on var.params
    basic.params <- csv2list(param.file(data.directory, file.id))
    var.params <- read.csv(var.params.file(data.directory, file.id))
    these.params <- modifyParams(basic.params, var.params, indices[["param"]])
    
    #' Generate a list of simulation labels unique for each parameter index (if not supplied)
    if (is.null(sim.labels)) sim.labels <- (indices[["node"]] - 1) * sims.per.node + seq_len(sims.per.node)
    #' Generate a unique seed for each simulation run
    seeds <- data.frame(sim.label = sim.labels, seed = sample(1:.Machine$integer.max, length(sim.labels)))
    
    #' Filename for the results and seeds
    output.filename <- output.file(data.directory, file.id, indices[["param"]], indices[["node"]])
    seed.filename <- seed.file(data.directory, file.id, indices[["param"]], indices[["node"]])
    
    #' Determine which simulation type to run (one of c("survival", ...))
    if (output == "survival") runSurvival(cell.type, these.params, seeds, output.filename, seed.filename, local.parallel)
    if (output == "growth") runGrowth(cell.type, these.params, seeds, output.filename, seed.filename, local.parallel)
    if (output == "molecules") runMolecules(cell.type, these.params, seeds, output.filename, seed.filename, local.parallel)
}


#' Compute lineage survival probabilities
runSurvival <- function(cell.type, these.params, seeds, output.filename, seed.filename, local.parallel) {
    #' Initialise our single cell that initiates the lineage
    initial.cell <- c(1)
    names(initial.cell) <- cell.type
    #' Length of the experiment
    t.max <- 1200
    #' Function to run to extract the population size at the end of each simulation
    getSurvival <- function(sim.label) {
        #' Set seed
        set.seed(seeds[seeds$sim.label == sim.label, "seed"])
        #' Generate an initial condition
        initial.cell.data <- simulateIC(these.params, initial.cell)
        #' Run the simulation over the duration of the experiment
        out <- simulateCells(these.params, initial.cell.data, cell.type, t.max, max.num.cells = 200, store.info = FALSE)
        #' Calculate population size curve
        pop.size <- analysePopSize(out$family.tree)
        #' Create a new output dataframe of the final size and exit time, using teh unique simulation label as an ID column
        pop.survive <- data.frame(
            sim.label = sim.label,
            num.cells = pop.size[nrow(pop.size), "Total"],
            exit.time = pop.size[nrow(pop.size), "time"]
        )
        #' Add survival column
        pop.survive$survive <- as.integer(pop.survive$num.cells > 0)
        return(pop.survive)
    }
    
    #' Run simulations within loop
    if (local.parallel) {
        pop.survival <- parallel.DFapply(seeds$sim.label, function(sim.label) {
            tryCatch(getSurvival(sim.label),
                     error = function(e) NULL)
        })
        #' Store final population sizes
        write.csv(pop.survival, file = output.filename, row.names = F)
        #' Store seeds
        write.csv(seeds[seeds$sim.label %in% pop.survival$sim.label, ], file = seed.filename, row.names = F)
    } else {
        for (sim.label in seeds$sim.label) {
            #' Run the simulation for each label
            pop.survival <- tryCatch(getSurvival(sim.label),
                                     error = function(e) NULL)
            #' Store final population sizes and seeds (appending to file each time)
            if (!is.null(pop.survival)) {
                write.table(pop.survival, output.filename, sep = ",", col.names = !file.exists(output.filename), append = T, row.names = F)
                write.table(data.frame(sim.label = sim.label, seed = seeds[seeds$sim.label == sim.label, "seed"]), seed.filename, sep = ",", col.names = !file.exists(seed.filename), append = T, row.names = F)
            }
        }
    }
}

#' Compute growth rates
runGrowth <- function(cell.type, these.params, seeds, output.filename, seed.filename, local.parallel) {
    #' Initialise 1000 cells per simulation
    initial.cells <- c(1000)
    names(initial.cells) <- cell.type
    #' Length of the experiment = 3hrs
    t.max <- 180
    #' Function to run to extract the population size during the simulations
    getGrowth <- function(sim.label) {
        #' Set seed
        set.seed(seeds[seeds$sim.label == sim.label, "seed"])
        #' Generate initial conditions
        initial.cell.data <- simulateIC(these.params, initial.cells)
        #' Run the simulation over the duration of the experiment
        out <- simulateCells(these.params, initial.cell.data, rep.int(cell.type, nrow(initial.cell.data)), t.max, max.num.cells = Inf, store.info = FALSE)
        #' Calculate population size curve, and simplify
        pop.size <- analysePopSize(out$family.tree)
        pop.size <- data.frame(
            sim.label = sim.label,
            sample.time = pop.size$time,
            num.cells = pop.size$Total
        )
        #' Now sample the population size at 10 minute intervals
        pop.size.sampled <- DFapply(seq(0, max(pop.size$sample.time), by = 10), function(t) {
            tmp <- pop.size[pop.size$sample.time <= t,]
            tmp <- tmp[tmp$sample.time == max(tmp$sample.time), ]
            tmp$sample.time <- t
            return(tmp)
        })
        return(pop.size.sampled)
    }
    
    #' Run simulations within loop
    if (local.parallel) {
        pop.size <- parallel.DFapply(seeds$sim.label, function(sim.label) {
            tryCatch(getGrowth(sim.label),
                     error = function(e) NULL)
        })
        #' Store final population sizes
        write.csv(pop.size, file = output.filename, row.names = F)
        #' Store seeds
        write.csv(seeds[seeds$sim.label %in% pop.size$sim.label, ], file = seed.filename, row.names = F)
    } else {
        for (sim.label in seeds$sim.label) {
            #' Run the simulation for each label
            pop.size <- #tryCatch(
                getGrowth(sim.label)#,
                        #         error = function(e) NULL)
            #' Store final population sizes and seeds (appending to file each time)
            if (!is.null(pop.size)) {
                write.table(pop.size, output.filename, sep = ",", col.names = !file.exists(output.filename), append = T, row.names = F)
                write.table(data.frame(sim.label = sim.label, seed = seeds[seeds$sim.label == sim.label, "seed"]), seed.filename, sep = ",", col.names = !file.exists(seed.filename), append = T, row.names = F)
            }
        }
    }
}


#' Compute molecules distributions
runMolecules <- function(cell.type, these.params, seeds, output.filename, seed.filename, local.parallel) {
    #' Initialise our single cell that initiates the lineage
    initial.cell <- c(1)
    names(initial.cell) <- cell.type
    #' Function to run to extract the number of molecules
    getMolecules <- function(sim.label) {
        #' Set seed
        set.seed(seeds[seeds$sim.label == sim.label, "seed"])
        #' Simulate a single cell for this duration using the initial condition generator
        cell.data <- simulateIC(these.params, initial.cell, t.max = 1200)
        #' Add a sim.label column
        cell.data$sim.label <- sim.label
        return(cell.data)
    }
    #' Run simulations within loop
    if (local.parallel) {
        molecules <- parallel.DFapply(seeds$sim.label, function(sim.label) {
            tryCatch(getMolecules(sim.label),
                     error = function(e) NULL)
        })
        #' Store final population sizes
        write.csv(molecules, file = output.filename, row.names = F)
        #' Store seeds
        write.csv(seeds[seeds$sim.label %in% molecules$sim.label, ], file = seed.filename, row.names = F)
    } else {
        for (sim.label in seeds$sim.label) {
            #' Run the simulation for each label
            molecules <- tryCatch(getMolecules(sim.label),
                                     error = function(e) NULL)
            #' Store final population sizes and seeds (appending to file each time)
            if (!is.null(molecules)) {
                write.table(molecules, output.filename, sep = ",", col.names = !file.exists(output.filename), append = T, row.names = F)
                write.table(data.frame(sim.label = sim.label, seed = seeds[seeds$sim.label == sim.label, "seed"]), seed.filename, sep = ",", col.names = !file.exists(seed.filename), append = T, row.names = F)
            }
        }
    }
}
