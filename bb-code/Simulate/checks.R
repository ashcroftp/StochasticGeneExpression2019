#' checks.R
#' Author: Peter Ashcroft, ETH Zurich

checkConfigFile <- function(config.file) {
    if (!file.exists(config.file)) {
        stop(paste0("config file does not exist: ", config.file))
    }
}

checkDrug <- function(drug) {
    if (length(drug) > 1) {
        stop(paste0("Only one drug at a time: ", paste0(drug, collapse = ", ")))
    }
    if (!drug %in% c("cefotaxime", "ciprofloxacin", "rifampicin")) {
        stop(paste0("Drug not permitted: ", drug))
    }
}

checkCellType <- function(cell.type) {
    if (length(cell.type) > 1) {
        stop(paste0("Only one cell type at a time: ", paste0(cell.type, collapse = ", ")))
    }
    if (!cell.type %in% c("KO", "WT", "REG-ON", "REG-OFF", "REG-BURST", "STRUCT-BIND", "STRUCT-CAT")) {
        stop(paste0("Cell type not permitted: ", cell.type))
    }
}

checkOutput <- function(output) {
    if (length(output) > 1) {
        stop(paste0("Only one output at a time: ", paste0(output, collapse = ", ")))
    }
    if (!output %in% c("survival", "growth", "molecules")) {
        stop(paste0("Output format not permitted: ", output))
    }
}

checkTotalSims <- function(total.sims) {
    if (length(total.sims) > 1) {
        stop(paste0("total.sims must be a single integer: ", paste0(total.sims, collapse = ", ")))
    }
    if (!total.sims > 0) {
        stop(paste0("The number of simulations must be greater than zero: ", total.sims))
    }
}

checkSimsPerNode <- function(sims.per.node, total.sims) {
    if (length(sims.per.node) > 1) {
        stop(paste0("sims.per.node must be a single integer: ", paste0(sims.per.node, collapse = ", ")))
    }
    if (sims.per.node < 0) {
        stop(paste0("The number of simulations per node cannot be less than zero: ", total.sims))
    }
    if (sims.per.node > total.sims) {
        stop(paste0("The number of simulations per node must be less than or equal to the total number of simulations: ", sims.per.node, " vs ", total.sims))
    }
}

checkVarParams <- function(var.params, allowed.names) {
    #' Make sure we modify parameters that exist
    allowed.names.MULT <- paste0(allowed.names, ".MULT")
    false.pars <- !sapply(names(var.params), function(param) param %in% c(allowed.names, allowed.names.MULT))
    if (any(false.pars)) {
        stop(paste0("Unknown parameters: ", paste0(names(var.params)[false.pars], collapse = ", ")))
    }
    #' Are the values of the modified parameters allowed?
    for (param in names(var.params)) {
        param.values <- unique(var.params[, param])
        if (any(param.values < 0) & !(param %in% c("rate.minGrowth", "rate.minGrowth.MULT", "fitness.cost", "fitness.cost.MULT"))) {
            stop(paste0(param, " cannot be less than zero: ", paste0(param.values, collapse = ", ")))
        }
        if (param == "MIC.fraction") {
            if (any(param.values <= 0) | any(param.values > 1)) {
                stop(paste0("mic.fraction must be greater than zero and less than or equal to one: ", paste0(param.values, collapse = ", ")))
            }
        }
    }
    rm(param)
}

checkFileID <- function(file.id) {
    if (length(checkFileID) > 1) {
        stop(paste0("file.id should be a single name: ", paste0(file.id, collapse = ", ")))
    }
    if (!nchar(file.id) > 0) {
        stop(paste0("file.id should contain at least one character: ", file.id))
    }
}

checkDataDirectory <- function(data.directory) {
    if (length(data.directory) > 1) {
        stop(paste0("data.directory should be a single name: ", paste0(data.directory, collapse = ", ")))
    }
    if (!dir.exists(data.directory)) {
        stop(paste0("data.directory does not exist: ", data.directory))
    }
}
