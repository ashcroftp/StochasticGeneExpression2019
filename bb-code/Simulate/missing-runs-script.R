#' missing-runs-script.R
#' Author: Peter Ashcroft, ETH Zurich

#' Determine which simulations terminated early, and re-run these
#' Dependencies: ----
source("Simulate/simulate-functions.R")

#' Get arguments from terminal
args <- commandArgs(trailingOnly = TRUE)

#' Script to run
if (length(args) == 0) {
    stop("You must pass a config file as the argument (no arguments passed)")
}

if (length(args) == 1) {
    config.file <- args[1]
    #' Check if the config file exists
    checkConfigFile(config.file)
    #' Source the config file
    source(config.file)
    #' Compute number of nodes used and number of parameters
    if (sims.per.node == 0) sims.per.node <- total.sims
    num.nodes <- total.sims/sims.per.node
    num.params <- nrow(var.params)
    
    #' First identify if any output files are missing
    missing.files <- getMissingFiles(data.directory, file.id, num.params, num.nodes)
    if (nrow(missing.files) > 0) {
        print("Missing files:")
        print(missing.files)
    } else {
        print("No files missing.")
        #' Now check to see if any sims are missing
        missing.sims <- getMissingSims(data.directory, file.id, num.params, num.nodes, sims.per.node)
        if (length(missing.sims) > 0) {
            print("Missing simulations:")
            print(sapply(missing.sims, function(k) k$sim.labels, simplify = F))
        } else {
            print("No simulations missing.")
        }
    }
    
    
} else if (length(args) == 2) {
    config.file <- args[1]
    #' Check if the config file exists
    checkConfigFile(config.file)
    #' Source the config file
    source(config.file)
    #' Compute number of nodes used and number of parameters
    if (sims.per.node == 0) sims.per.node <- total.sims
    num.nodes <- total.sims/sims.per.node
    num.params <- nrow(var.params)
    #' Make sure second argument is run
    if (args[2] == "run") {
        #' First we identify if files are missing, and then rerun these sims
        missing.files <- getMissingFiles(data.directory, file.id, num.params, num.nodes)
        if (nrow(missing.files) > 0) {
            print(paste0("A total of ", nrow(missing.files), " files are missing." ))
            print("Now resubmitting these missing simulations.")
            #' Re-submit these simulations
            invisible(
                lapply(missing.files$index, function(index) {
                    system(
                        paste(
                            "bsub", "-W", queue,
                            "\"R --vanilla --slave < Simulate/simulate-script.R --args",
                            data.directory, file.id, cell.type, output, num.params, num.nodes, sims.per.node, as.integer(index), "\"",
                            sep = " "
                        )
                    )
                })
            )
        } else {
            print("No files missing, now checking for missing simulations")
            #' Identify which simulation runs are missing
            missing.sims <- getMissingSims(data.directory, file.id, num.params, num.nodes, sims.per.node)
            #' Re-submit these simulations
            invisible(
                lapply(names(missing.sims), function(index) {
                    system(
                        paste(
                            "bsub", "-W", queue,
                            "\"R --vanilla --slave < Simulate/simulate-script.R --args",
                            data.directory, file.id, cell.type, output, num.params, num.nodes, sims.per.node, as.integer(index), paste(missing.sims[[index]]$sim.labels, collapse = " "), "\"",
                            sep = " "
                        )
                    )
                })
            )
        }
        
    } else {
        stop(paste0("The second argument is invalid: ", args[2]))
    }
} else {
    stop(paste0("Incorrect number of arguments passed: ", length(args)))
}
