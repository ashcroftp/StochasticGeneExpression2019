#' simulate-population.R
#' Author: Peter Ashcroft, ETH Zurich

#' Here we store the functions that perform simulations of populations of a single cell type.
#' The name of a configuration file must be passed, and then simulations based on this file will be run (on EULER or not).

#' Working directory to base (as we source multiple functions, that rely on being in that directory)
wd <- basename(getwd())
while (wd != "bb-code") {
    setwd("../")
    wd <- basename(getwd())
}

#' Dependencies: ----
source("R/basicFunctions.R")
source("R/simulationFunctions.R")
source("R/dataAnalysisFunctions.R")
source("Simulate/checks.R")
source("Simulate/filenames.R")
#' Enable byte-code compiling for a little speed up
library(compiler)
enableJIT(1)

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
runSim <- function(cell.type, output, num.realisations, data.directory, file.id, var.param.index, node.index, local.parallel = FALSE) {
    #' Load parameters from local file, and modify parameters based on var.params
    basic.params <- csv2list(param.file(data.directory, file.id))
    var.params <- read.csv(var.params.file(data.directory, file.id))
    these.params <- modifyParams(basic.params, var.params, var.param.index)
    
    if (output == "survival") {
        #' Calculate survival probabilities
        #' 
        #' Assign and save seeds for each simulation
        seeds <- sample(1:.Machine$integer.max, num.realisations)
        list2csv(as.list(seeds), seed.file(data.directory, file.id, var.param.index, node.index))
        #' Initialise our single cell
        initial.cell <- c(1)
        names(initial.cell) <- cell.type
        #' Length of the experiment
        t.max <- 1200
        #' Function to run to extract the population size at the end of each simulation
        runSurvival <- function(run.label) {
            #' Set seed
            set.seed(seeds[run.label])
            #' Initialise a single WT cell
            initial.cell.data <- simulateIC(these.params, initial.cell)
            #' Run
            out <- simulateCells(these.params, initial.cell.data, cell.type, t.max, max.num.cells = 100, store.info = FALSE)
            #' Calculate population size curve
            pop.size <- analysePopSize(out$family.tree)
            #' Create a new output dataframe with a unique simulation label
            pop.survive <- data.frame(
                sim.label = max(node.index - 1, 0) * num.realisations + run.label,
                num.cells = pop.size[nrow(pop.size), "Total"],
                exit.time = pop.size[nrow(pop.size), "time"]
            )
            #' Add survival column
            pop.survive$survive <- as.integer(pop.survive$num.cells > 0)
            return(pop.survive)
        }
        #' Run simulations within loop
        if (local.parallel) {
            pop.survival <- parallel.DFapply(seq_len(num.realisations), runSurvival)
        } else {
            pop.survival <- DFapply(seq_len(num.realisations), runSurvival)
        }
        #' Store final population sizes
        write.csv(pop.survival, file = output.file(data.directory, file.id, var.param.index, node.index), row.names = FALSE)
        
    } else if (output == "growth") {
        #' Output popluation growth curves of 100 initial cells over the course of 1 hour
        #' 
        #' Assign and save seeds for each simulation
        seeds <- sample(1:.Machine$integer.max, num.realisations)
        list2csv(as.list(seeds), seed.file(data.directory, file.id, var.param.index, node.index))
        #' Initialise our cells
        initial.cells <- c(1000)
        names(initial.cells) <- cell.type
        #' Length of the experiment
        t.max <- 60
        #' Function to run to extract the population size at the end of each simulation
        runGrowth <- function(run.label) {
            #' Set seed
            set.seed(seeds[run.label])
            #' Initialise a single WT cell
            initial.cell.data <- simulateIC(these.params, initial.cells)
            #' Run
            out <- simulateCells(these.params, initial.cell.data, rep.int(cell.type, nrow(initial.cell.data)), t.max, max.num.cells = Inf, store.info = FALSE)
            #' Calculate population size curve
            pop.size <- analysePopSize(out$family.tree)
            #' Return cleaner dataframe wuth unique simulation label
            data.frame(
                sim.label = max(node.index - 1, 0) * num.realisations + run.label,
                sample.time = pop.size$time,
                num.cells = pop.size$Total
            )
        }
        #' Run simulations within loop
        if (local.parallel) {
            growth.curves <- parallel.DFapply(seq_len(num.realisations), runGrowth)
        } else {
            growth.curves <- DFapply(seq_len(num.realisations), runGrowth)
        }
        #' Store final population sizes
        write.csv(growth.curves, file = output.file(data.directory, file.id, var.param.index, node.index), row.names = FALSE)
    } else if (output == "molecules") {
        #' Return the number of each within-cell molecule for a number of simulations
        #' 
        #' Assign and save seed for this simulation
        seed <- sample(1:.Machine$integer.max, 1)
        list2csv(as.list(seed), seed.file(data.directory, file.id, var.param.index, node.index))
        #' Initialise our cells
        cells <- c(num.realisations)
        names(cells) <- cell.type
        #' Simulate multiple cells for multiple generations, and just retain the last timepoint
        #' Here we can use the initial-condition function, which does just this
        set.seed(seed)
        cell.data <- simulateIC(these.params, cells, t.max = 1200)
        #' Store data
        write.csv(cell.data, file = output.file(data.directory, file.id, var.param.index, node.index), row.names = FALSE)
    }
}



#' Read in command-line arguments and run the simulations.
#' If only the config file is passed, then parameters are checked and the simulations are triggered.
#' If config file and another argument is passed, then we execute the second argument
#' If multiple arguments are passed, then the simulations are executed, either locally or on EULER.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
    stop("You must pass a config file as the argument (no arguments passed)")
}
if (length(args) == 1) {
    config.file <- args[1]
    #' Check if the config file exists
    checkConfigFile(config.file)
    #' Source the config file
    source(config.file)
    #' Check if any of the required parameters are missing
    required.pars <- c("drug", "cell.type", "output", "total.sims", "sims.per.node", "var.params", "file.id", "data.directory")
    missing.pars <- !sapply(required.pars, exists)
    if (any(missing.pars)) {
        stop(paste0("Missing parameters: ", paste0(required.pars[missing.pars], collapse = ", ")))
    }
    #' Check the compatibility of the parameters
    checkDrug(drug)
    checkCellType(cell.type)
    checkOutput(output)
    checkTotalSims(total.sims)
    checkSimsPerNode(sims.per.node, total.sims)
    #' Do we use the cluster?
    EULER <- ifelse(sims.per.node == 0, FALSE, TRUE)
    #' Check compatibility of files and existence of directories
    checkFileID(file.id)
    checkDataDirectory(data.directory)
    #' Load the parameters, and check the var.param dataframe is compatible
    basic.params <- loadParameters(drug)
    checkVarParams(var.params, names(basic.params))
    
    #' Locally save the parameter list
    list2csv(basic.params, param.file(data.directory, file.id))
    #' Save the parameter combinations over which we will simulate
    write.csv(var.params, file = var.params.file(data.directory, file.id), row.names = F)
    #' Save the working environment
    sink(file = info.file(data.directory, file.id))
    print(
        list(
            drug = drug,
            cell.type = cell.type,
            output = output,
            total.sims = total.sims,
            sims.per.node = sims.per.node,
            EULER = EULER,
            file.id = file.id,
            data.directory = data.directory,
            timestamp = timestamp(quiet = T),
            version = R.Version()
        )
    )
    sink()
    
    if (EULER) {
        #' Run simulations on EULER using the following number of nodes
        num.nodes <- total.sims / sims.per.node
        invisible(
            lapply(seq_len(nrow(var.params)), function(parameter.index) {
                system(
                    paste(
                        paste("bsub -J \"", drug, cell.type, output, parameter.index, paste("[1-", num.nodes, "]\"", sep = ""), sep = " "),
                        paste("\"R --vanilla --slave < Simulate/simulate.R --args", cell.type, output, sims.per.node, data.directory, file.id, parameter.index, "\\$LSB_JOBINDEX\"", sep = " "),
                        sep = " "
                    )
                )
            })
        )
    }
    else{
        invisible(
            lapply(seq_len(nrow(var.params)), function(parameter.index) {
                runSim(cell.type = cell.type, output = output, num.realisations = total.sims, data.directory = data.directory, file.id = file.id, var.param.index = parameter.index, node.index = 0, local.parallel = TRUE)
            })
        )
    }
} else if (length(args) == 2) {
    config.file <- args[1]
    #' Check if the config file exists
    checkConfigFile(config.file)
    #' Source the config file
    source(config.file)
    #' Make sure second argument is compress
    if (args[2] == "compress") {
        #' Parameter indices
        parameter.indices <- seq_len(nrow(var.params))
        #' Loop over parameter indices, load the data, save the good stuff
        data <- DFapply(parameter.indices, function(parameter.index) {
            #' Identify the output files
            output.files <- Sys.glob(output.file(data.directory, file.id, parameter.index, "*"))
            seeds.files <- Sys.glob(seed.file(data.directory, file.id, parameter.index, "*"))
            #' Now load each file, and combine into single DF
            data <- DFapply(seq_along(output.files), function(node.index) {
                output.data <- read.csv(output.files[node.index])
                seeds <- csv2list(seeds.files[node.index])
                #' Add the metadata to the output
                #output.data$seed <- unlist(seeds)
                output.data$parameter.index <- parameter.index
                for (i in names(var.params)) {
                    output.data[,i] <- var.params[parameter.index, i]
                }
                rm(i)
                return(output.data)
            })
            return(data)
        })
        #' Now write to file
        write.csv(data, file = combined.file(data.directory, file.id), row.names = F)
        #' Move raw data to a new directory, and then archive it
        raw.dir <- paste0(data.directory, file.id)
        dir.create(raw.dir)
        system(paste("mv", output.file(data.directory, file.id, "*", "*"), raw.dir, sep = " "))
        system(paste("mv", seed.file(data.directory, file.id, "*", "*"), raw.dir, sep = " "))
        system(paste("tar -czf", paste0(raw.dir, ".tar.gz"), "-C", raw.dir, ".", sep = " "))
    } else {
        stop(paste0("The second argument is invalid: ", args[2]))
    }
} else if (length(args) == 7) {#' We run the simulation on EULER
    #' Cast arguments as integer, numeric or character
    cell.type <- args[1]
    output <- args[2]
    num.realisations <- as.integer(args[3])
    data.directory <- args[4]
    file.id <- args[5]
    var.param.index <- as.integer(args[6])
    node.index <- as.integer(args[7])
    #' Run the simulation
    runSim(cell.type, output, num.realisations, data.directory, file.id, var.param.index, node.index)
} else {
    stop(paste0("Incorrect number of arguments passed: ", length(args)))
}
