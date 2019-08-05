#' compress-script.R
#' Author: Peter Ashcroft, ETH Zurich

#' This script is designed to be run from the command line using the syntax 
#' \code{R --vanilla --slave < Simulate/compress-script.R --args config.R}
#' The name of a configuration file must be passed, and then the output files
#' from this will be combined into a single output file.

#' Dependencies: ----
source("Simulate/simulate-functions.R")

#' Read in command-line arguments and compress the output data
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
    #' Number of parameter combinations
    num.params <- nrow(var.params)
    #' First we loop over the parameter indices to compress the files from each node into a single file
    for (param.index in seq_len(num.params)) {
        #' Identify the output and seed files for this index
        output.files <- Sys.glob(output.file(data.directory, file.id, param.index, "*"))
        seed.files <- Sys.glob(seed.file(data.directory, file.id, param.index, "*"))
        
        #' Now load each output file, combine into single DF, and then save as a single file
        data <- DFapply(output.files, function(output.filename) {
            output.data <- read.csv(output.filename)
            #' Add the parameter index to DF
            output.data$param.index <- param.index
            #' Add the var.params to the DF
            for (i in names(var.params)) {
                output.data[,i] <- var.params[param.index, i]
            }
            rm(i)
            return(output.data)
        })
        write.csv(data, file = output.file(data.directory, file.id, param.index), row.names = F)
        
        #' Repeat for seeds
        seeds <- DFapply(seed.files, function(seed.filename) {
            seed.data <- read.csv(seed.filename)
            #' Add the parameter index to DF
            seed.data$param.index <- param.index
            return(seed.data)
        })
        write.csv(seeds, file = seed.file(data.directory, file.id, param.index), row.names = F)
        
        #' Remove old files
        for (output.filename in output.files) system(paste("rm", output.filename, sep = " "))
        for (seed.filename in seed.files) system(paste("rm", seed.filename, sep = " "))
        rm(output.filename, seed.filename)
    }
    rm(param.index)
    
    #' Now we combine these outputs using terminal commands, rather than reading in the data
    combineFiles <- function(data.directory, result.files, output.filename, chunk = 1000) {
        #' Divide the filenames into chunks of size \code{chunk}
        result.files <- split(result.files, ceiling(seq_along(result.files)/chunk))
        #' Create a blank filename
        system(paste("touch", output.filename, sep = " "))
        for (i in names(result.files)) {
            #' Concatenate the chunk of results files with the output file, and temporariliy store
            system(paste("awk 'FNR>1 || NR==1'", output.filename, paste(result.files[[i]], collapse = " "), ">", paste0(data.directory, "temp.dat"), "; echo", i , sep = " "))
            #' Rename the temp file to the output file
            system(paste("mv", paste0(data.directory, "temp.dat"), output.filename, sep = " "))
            #' Remove the files that have been concatenated
            system(paste("rm", paste(result.files[[i]], collapse = " "), sep = " "))
        }
    }
    #' Run the combining function for outputs and seeds
    combineFiles(data.directory, sapply(seq_len(num.params), function(param.index) output.file(data.directory, file.id, param.index)), output.file(data.directory, file.id))
    combineFiles(data.directory, sapply(seq_len(num.params), function(param.index) seed.file(data.directory, file.id, param.index)), seed.file(data.directory, file.id))
} else {
    stop(paste0("Too many arguments: ", length(args)))
}
