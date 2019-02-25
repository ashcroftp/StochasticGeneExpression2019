#' filenames.R
#' Author: Peter Ashcroft, ETH Zurich

#' Defines the filenames globally

#' How are numbers added to filenames
num.format <- function(n) ifelse(is.numeric(n), sprintf("%03d", n), n)

#' Filename functions
param.file <- function(data.directory, file.id) paste0(data.directory, file.id, ".params")
var.params.file <- function(data.directory, file.id) paste0(data.directory, file.id, ".varparams")
info.file <- function(data.directory, file.id) paste0(data.directory, file.id, ".info")
seed.file <- function(data.directory, file.id, parameter.index, node.index) paste0(data.directory, file.id, "_seeds_", num.format(parameter.index), "_", num.format(node.index), ".dat")
output.file <- function(data.directory, file.id, parameter.index, node.index) paste0(data.directory, file.id, "_output_", num.format(parameter.index), "_", num.format(node.index), ".dat")
combined.file <- function(data.directory, file.id) paste0(data.directory, file.id, "_combined", ".dat")