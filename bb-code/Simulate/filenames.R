#' filenames.R
#' Author: Peter Ashcroft, ETH Zurich

#' Defines the filenames globally

#' How are numbers added to filenames
num.format <- function(n) ifelse(is.numeric(n), sprintf("%05d", n), n)

## Filename functions ----
#' Basic model parameters
param.file <- function(data.directory, file.id) {
    paste0(data.directory, file.id, ".params")
}

#' Parameters we vary to screen over (i.e. drug concentration)
var.params.file <- function(data.directory, file.id) {
    paste0(data.directory, file.id, ".varparams")
}

#' System state at time of simulations
info.file <- function(data.directory, file.id) {
    paste0(data.directory, file.id, ".info")
}

#' Seeds of each simulation run (with or without node index following compression)
seed.file <- function(data.directory, file.id, param.index = NULL, node.index = NULL) {
    if (is.null(param.index) & is.null(node.index)) {
        paste0(data.directory, file.id, "_seeds", ".dat")
    } else if (is.null(node.index)) {
        paste0(data.directory, file.id, "_seeds_", num.format(param.index), ".dat")
    } else {
        paste0(data.directory, file.id, "_seeds_", num.format(param.index), "_", num.format(node.index), ".dat")
    }
}

#' Data from each simulation run (with or without node index following compression)
output.file <- function(data.directory, file.id, param.index = NULL, node.index = NULL) {
    if (is.null(param.index) & is.null(node.index)) {
        paste0(data.directory, file.id, "_output", ".dat")
    } else if (is.null(node.index)) {
        paste0(data.directory, file.id, "_output_", num.format(param.index), ".dat")
    } else {
        paste0(data.directory, file.id, "_output_", num.format(param.index), "_", num.format(node.index), ".dat")
    }
}

#' Combined data file across all parameters following compression
combined.file <- function(data.directory, file.id) paste0(data.directory, file.id, "_combined", ".dat")
