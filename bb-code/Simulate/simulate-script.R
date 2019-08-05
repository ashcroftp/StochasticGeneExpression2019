#' simulate-script.R
#' Author: Peter Ashcroft, ETH Zurich

#' This script is designed to be run from the command line using the syntax 
#' \code{R --vanilla --slave < Simulate/simulate-script.R --args config.R}
#' The name of a configuration file must be passed, and then simulations based on this file will be run (on EULER or not).

#' Dependencies: ----
source("Simulate/simulate-functions.R")
#' Enable byte-code compiling for a little speed up
library(compiler)
enableJIT(1)

#' Read in command-line arguments and run the simulations.
#' If only the config file is passed, then parameters are checked and the simulations are triggered.
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
    #' If we use EULER, have we specified a correct queue
    if(EULER) checkQueue(queue)
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
            queue = queue,
            file.id = file.id,
            data.directory = data.directory,
            timestamp = timestamp(quiet = T),
            version = R.Version()
        )
    )
    sink()
    
    #' Now we run the simulation functions either on EULER or locally
    if (EULER) {
        #' Run simulations on EULER using the following number of parameters and nodes, and hence the computed number of combinations
        num.params <- nrow(var.params)
        num.nodes <- total.sims / sims.per.node
        num.combinations <- num.params * num.nodes
        #' RUN
        #invisible(
        #    system(
        #        paste(
        #            paste("bsub -W 24:00 -J \"", file.id, paste("[1-", num.combinations, "]\"", sep = ""), sep = ""),
        #            paste("\"R --vanilla --slave < Simulate/simulate-script.R --args", data.directory, file.id, cell.type, output, num.params, num.nodes, sims.per.node, "\\$LSB_JOBINDEX\"", sep = " "),
        #            sep = " "
        #        )
        #    )
        #)
        
        invisible(
            system(
                paste(
                    "bsub", "-W" , queue, "-J",
                    paste("\"", file.id, "[1-", num.combinations, "]\"", sep = ""),
                    "\"R --vanilla --slave < Simulate/simulate-script.R --args", data.directory, file.id, cell.type, output, num.params, num.nodes, sims.per.node, "\\$LSB_JOBINDEX\"",
                    sep = " "
                )
            )
        )
    }
    else{
        #' Run everything locally
        num.params <- nrow(var.params)
        sims.per.node <- total.sims
        #' Run
        invisible(
            lapply(seq_len(num.params), function(index) {
                runSim(data.directory, file.id, cell.type, output, num.params, 1, sims.per.node, index, local.parallel = F)
            })
        )
    }
} else if (length(args) == 8) {
    #' Run the simulations, casting arguments as integer or character
    runSim(data.directory = args[1], file.id = args[2], cell.type = args[3], output = args[4], num.params = as.integer(args[5]),
           num.nodes = as.integer(args[6]), sims.per.node = as.integer(args[7]), index = as.integer(args[8]))
    
} else if (length(args) > 8) {
    #' Now we are supplied a list of simulation labels to repeat the simulations for, which we pass as a vector to runSim()
    #' Run the simulations, casting arguments as integer or character
    runSim(data.directory = args[1], file.id = args[2], cell.type = args[3], output = args[4], num.params = as.integer(args[5]),
           num.nodes = as.integer(args[6]), sims.per.node = as.integer(args[7]), index = as.integer(args[8]), sim.labels = as.integer(args[9:length(args)]))
} else {
    stop(paste0("Incorrect number of arguments passed: ", length(args)))
}
