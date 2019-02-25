#' basicFunctions.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' Some basic functions that I will repeatedly use,
#' so I will define them here

#' Dependencies:
library(reshape2)
library(parallel)

#' Function list:
#' - parallel.lapply()
#' - DFapply()
#' - parallel.DFapply()
#' - splitFilenames()
#' - list2csv()
#' - csv2list()


#' A quick and easy way of parallelising lapply
#' @param X (vector or list): What is to be looped over.
#' @param FUN (function): First argument must be the value of \code{X}.
#' @param ...: Further arguments to function
#' @return List of outputs from \code{FUN} applied to each element of X.
parallel.lapply <- function(X, FUN, ...){
    #' First initialise a cluster
    #' We need to know how many cores we have, and set the number we can use
    #' to one less than this
    num.cores <- max(1, detectCores() - 1)
    #' Initialise the cluster
    cl <- makeCluster(num.cores, type = "FORK")
    #' Set a seed for the random number generator, which then seeds all new jobs
    #' to ensure no random sequences are used on multiple cores
    clusterSetRNGStream(cl, 111)
    #' Now our function
    output <- parLapply(cl, X, FUN, ...)
    #' Finally, close the cluster we created
    stopCluster(cl)
    return(output)
}

#' Sometimes we want to use lapply to bind dataframes,
#' so I define the function \code{DFapply} to do this for me
#' @param X (vector or list): What is to be looped over.
#' E.g. the list of dataframes we want to bind
#' @param FUN (function): First argument must be the value of \code{X}. Must return a dataframe.
#' @param ...: Further arguments to function
#' @return Dataframe of the combined dataframes from \code{FUN}.
DFapply <- function(X, FUN, ...) {
    df.list <- lapply(X, FUN, ...)
    df <- do.call(rbind, df.list)
    return(df)
}

#' We can also parallelise the \code{DFapply} function.
#' @param X (vector or list): What is to be looped over.
#' E.g. the list of dataframes we want to bind
#' @param FUN (function): First argument must be the value of \code{X}. Must return a dataframe.
#' @param ...: Further arguments to function
#' @return Dataframe of the combined dataframes from \code{FUN}.
parallel.DFapply <- function(X, FUN, ...) {
    df.list <- parallel.lapply(X, FUN, ...)
    df <- do.call(rbind, df.list)
    return(df)
}


#' Split a file name into its constituent parts
#' @param filenames (character vector): List of filenames
#' @param col.names (character vector): What does each element of the filename correspond too.
#' @param split (character): What character do we split the file name about.
#' @return Dataframe of split filenames, with names as given by \code{col.names}
splitFilenames <- function(filenames, col.names = NULL, split = "_") {
    # Separate the file extension
    labels <- do.call(rbind, strsplit(filenames, split = ".", fixed = TRUE))
    if (ncol(labels) > 2) {
        stop("Too many dots (.) in the filenames.")
    } else {
        labels <- labels[,1]
    }
    # Split remaining parts
    labels <- data.frame(do.call(rbind, strsplit(labels, split = split, fixed = TRUE)))
    # Add names
    if (!is.null(col.names)) {
        if (ncol(labels) != length(col.names) ) {
            stop("Incorrect number of names in col.names")
        } else {
            names(labels) <- col.names
        }
    }
    return(labels)
}

#' Save a list as a csv file.
#' @param x (list): A named, 1D list. If no names, we assign index as the name
#' @param filename (character): The name of the file to save to
list2csv <- function(x, filename) {
    if (is.null(names(x))) names(x) <- seq_along(x)
    df <- DFapply(names(x), function(name) data.frame(name, x[[name]]))
    #write.csv(df, file = filename, row.names = F)
    write.table(df, sep = ",", file = filename, row.names = F, col.names = F)
}

#' Load a csv file into a list.
#' @param filename (character): The name of the file to load
csv2list <- function(filename) {
    df <- read.csv(filename, colClasses = c("character", "numeric"), header = F)
    x <- as.list(df[, 2])
    names(x) <- df[, 1]
    return(x)
}
