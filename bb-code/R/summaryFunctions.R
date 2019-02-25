#' summaryFunctions.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:

#' Function list:
#' - timeWeightedMean()
#' - timeWeightedMean.allVars()
#' - startVals.allVars()
#' - endVals.allVars()
#' - extractCellSummaryStatistics()
#' - extractCellSummaryStatistics.allCells()
#' - mean.by.generation()

## Time-weighted mean values per generation ----
timeWeightedMean <- function(timeseries){
    colnames(timeseries) <- c("time","x")
    n <- nrow(timeseries)
    
    time.start <- timeseries[1, "time"]
    time.end <- timeseries[n, "time"]
    
    time.diff <- data.frame(time=(timeseries[2:n, "time"]-timeseries[1:(n-1), "time"]), x=timeseries[1:(n-1), "x"])
    this.weighted.mean <- sum( (time.diff$time * time.diff$x) / (time.end - time.start) )
 
    return(this.weighted.mean)
}

timeWeightedMean.allVars <- function(cell.data){
    vars <- names(cell.data)[!names(cell.data) %in% c("time","alive","generation")]
    means <- sapply(vars, function(this.var) timeWeightedMean(cell.data[, c("time", this.var)]))
    names(means) <- paste("mean", names(means), sep=".")
    return(means)
}

## Values at the start and end of each generation ----
startVals.allVars <- function(cell.data){
    vars <- names(cell.data)[!names(cell.data) %in% c("time","alive","generation")]
    start.vals <- sapply(vars, function(this.var) cell.data[1, this.var] )
    names(start.vals) <- paste("start.val", names(start.vals), sep=".")
    return(start.vals)
}

endVals.allVars <- function(cell.data){
    vars <- names(cell.data)[!names(cell.data) %in% c("time","alive","generation")]
    end.vals <- sapply(vars, function(this.var) cell.data[nrow(cell.data), this.var] )
    names(end.vals) <- paste("end.val", names(end.vals), sep=".")
    return(end.vals)
}

extractCellSummaryStatistics <- function(cell.data){
    # Extract mean values for each variable
    mean.vals <- timeWeightedMean.allVars(cell.data)
    mean.vals <- data.frame(as.list(mean.vals) )
    
    # Extract start and end values for each variable (just after division and just before division)
    start.vals <- startVals.allVars(cell.data)
    start.vals <- data.frame(as.list(start.vals) )
    end.vals <- endVals.allVars(cell.data)
    end.vals <- data.frame(as.list(end.vals) )
    
    # Combine into single summary statistics 
    summary.statistics <- cbind(mean.vals, start.vals, end.vals)
    #summary.statistics <- cbind(start.vals, end.vals)
    
    # Add the generation label
    summary.statistics$generation <- cell.data[1, "generation"]
    
    return(summary.statistics)
}

extractCellSummaryStatistics.allCells <- function(model.params, cell.data.list){
    summary.statistics <- data.frame()
    for(i in 1:length(cell.data.list) ){
        thisCell.data <- cell.data.list[[as.character(i)]]
        
        # Add column for drug concentration
        cell.volume <- pi*(model.params$cell.diameter^2)*(2^(1+((thisCell.data[,"time"]-thisCell.data[1,"time"])/model.params$generationTime)))/4
        thisCell.data$Concentration <- thisCell.data$Drug / cell.volume
        
        # Compute summary statistics
        thisCell.summaryStatistics <- extractCellSummaryStatistics(thisCell.data)
        summary.statistics <- rbind(summary.statistics, thisCell.summaryStatistics)
    }
    
    return(summary.statistics)
}

means.by.generation <- function(summary.statistics){
    # Split data by generations (first generation cells in one dataframe, second generation cells...)
    summary.by.generation <- split(summary.statistics, summary.statistics$generation)
    # Take mean over all columns within each generation dataframe
    generation.means <- lapply(summary.by.generation, function(generation.stats) data.frame(as.list(apply(generation.stats, 2, mean) ) ) )
    #generation.means <- lapply(names(generation.means), function(this.generation) cbind(generation.means[[this.generation]], data.frame(generation=as.numeric(this.generation) ) ) )
    
    # Combine to single dataframe
    generation.means <- do.call("rbind", generation.means)
    
    return(generation.means)
}
