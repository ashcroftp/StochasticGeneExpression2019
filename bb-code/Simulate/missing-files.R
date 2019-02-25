#' missing-files.R
#' Author: Peter Ashcroft, ETH Zurich

source("Simulate/filenames.R")
source("Simulate/checks.R")

args <- commandArgs(trailingOnly = TRUE)
# <- c("config-rif-kappa.R")


getMissing <- function() {
    #' List the stored output files
    output.files <- paste0(data.directory, list.files(data.directory, pattern = "*_output_*"))
    #' List the expected output files
    expected <- expand.grid(node.index = seq_len(total.sims/sims.per.node), parameter.index = seq_len(nrow(var.params)) )
    expected <- expected[,c("parameter.index", "node.index")]
    expected$file <- sapply(seq_len(nrow(expected)), function(i) output.file(data.directory, file.id, expected[i, "parameter.index"], expected[i, "node.index"]))
    #' Which files are missing
    expected$missing <- sapply(seq_len(nrow(expected)), function(i) !expected[i, "file"] %in% output.files)
    return(expected[expected$missing, c("parameter.index", "node.index")])
}




if (length(args) == 0) {
    stop("You must pass a config file as the argument (no arguments passed)")
}
if (length(args) == 1) {
    config.file <- args[1]
    #' Check if the config file exists
    checkConfigFile(config.file)
    #' Source the config file
    source(config.file)
    #' Identify the missing files
    out <- getMissing()
    print(out)
} else if (length(args) == 2) {
    config.file <- args[1]
    #' Check if the config file exists
    checkConfigFile(config.file)
    #' Source the config file
    source(config.file)
    #' Make sure second argument is compress
    if (args[2] == "run") {
        #' Identify missing files
        missing <- getMissing()
        #' Re-run these simulations
        invisible(
            lapply(seq_len(nrow(missing)), function(i) {
                system(
                    paste("bsub", "\"R --vanilla --slave < Simulate/simulate.R --args", cell.type, output, sims.per.node, data.directory, file.id, missing[i, "parameter.index"], missing[i, "node.index"], "\"", sep = " ")
                )
            })
        )
        
    } else {
        stop(paste0("The second argument is invalid: ", args[2]))
    }
} else {
    stop(paste0("Incorrect number of arguments passed: ", length(args)))
}
