#' benchmark.R
#' Author: Peter Ashcroft, ETH Zurich

#' A script to time how long some simulations take

#' Working directory to base (as we source multiple functions, that rely on being in that directory)
wd <- basename(getwd())
while (wd != "bb-code") {
    setwd("../")
    wd <- basename(getwd())
}

#' Dependencies:
library(tictoc)

args <- commandArgs(trailingOnly = TRUE)

tic("simulation time")
system(paste("R --vanilla --slave < Simulate/simulate.R --args", args[1], sep = " "))
toc(log = T)
