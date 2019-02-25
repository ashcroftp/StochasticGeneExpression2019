#' meanFieldModel.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' ODE models of the within-cell processes.
#' This process is modelled in two ways.
#' First we continuoulsy dilute the contents of the cell across multiple generations.
#' Alternatively, we can integrate the ODES until the end of a generation, and then divide the contents in two for the daughter.
#' We call these 'dilute' and 'division' equations. 

#' Dependencies:
library(deSolve)

#' Function list:
#' - diluteEqns()
#' - integrateDiluteEqns()
#' - divisionEqns()
#' - integrateDivisionEqns()

#' Within-cell equations, where division is replaced by a continuous dilution
diluteEqns <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
        # Average area & volume of a cell
        SurfaceArea <- 2 * pi * (cell.diameter / 2) * cell.initialLength / log(2) + 2 * pi * (cell.diameter / 2) ^ 2
        Volume <- pi * (cell.diameter / 2) ^ 2 * cell.initialLength / log(2)
        # Drug concentration
        DrugConcentration <- drug / Volume
        
        # Quantities that are repeatedly used
        TargetBinding <- rate.target.binding * target.protein * DrugConcentration
        TargetUnbinding <- rate.dissociation * target.dimer
        TargetDeacylation <- rate.deacylation * target.dimer
        EffluxBinding <- rate.efflux.binding * efflux.protein * DrugConcentration
        EffluxUnbinding <- rate.dissociation * efflux.dimer
        EffluxDeacylation <- rate.deacylation * efflux.dimer
        Catalysis <- rate.catalysis * efflux.dimer
        DilutionRate <- log(2) / mean.generationTime
        
        # RHS of differential equations    
        d.target.protein    <- rate.target.protein.translation - TargetBinding + TargetUnbinding + TargetDeacylation - DilutionRate * target.protein
        d.efflux.DNA        <- rate.efflux.DNA.switchOn * (1 - efflux.DNA) - rate.efflux.DNA.switchOff * efflux.DNA
        d.efflux.mRNA       <- rate.efflux.mRNA.transcription * efflux.DNA - rate.efflux.mRNA.deg * efflux.mRNA - DilutionRate * efflux.mRNA
        d.efflux.protein    <- rate.efflux.protein.translation * efflux.mRNA - EffluxBinding + EffluxUnbinding + EffluxDeacylation + Catalysis - DilutionRate * efflux.protein
        d.drug              <- rate.drug.diffusion * SurfaceArea * (ifelse(t < drug.off, drug.outsideConcentration, 0) - DrugConcentration) - (TargetBinding + EffluxBinding) + (TargetUnbinding + EffluxUnbinding) - DilutionRate * drug
        d.target.dimer      <- TargetBinding - TargetUnbinding - TargetDeacylation - DilutionRate * target.dimer
        d.efflux.dimer      <- EffluxBinding - (EffluxUnbinding + EffluxDeacylation + Catalysis) - DilutionRate * efflux.dimer
        
        return(list(c(d.target.protein, d.efflux.DNA, d.efflux.mRNA, d.efflux.protein, d.drug, d.target.dimer, d.efflux.dimer)))
    })
}

integrateDiluteEqns <- function(t.max, params, IC){
    # Initial values for each variable
    vars <- c("target.protein", "efflux.DNA", "efflux.mRNA", "efflux.protein", "drug", "target.dimer", "efflux.dimer")
    init.cond <- unlist(IC[1, vars])
    # Times
    times <- seq.int(0, t.max)
    # Solve
    sol <- as.data.frame(ode(init.cond, times, diluteEqns, params))
    # Add volume
    sol$volume <- with(params, pi * (cell.diameter / 2) ^ 2 * cell.initialLength / log(2) )
    # Return
    return(sol)
}


#' Within-cell equations, where division is computed explicitly
divisionEqns <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
        # Average area & volume of a cell
        SurfaceArea <- 2 * pi * (cell.diameter / 2) * cell.initialLength * 2^((t %% mean.generationTime) / mean.generationTime) + 2 * pi * (cell.diameter / 2) ^ 2
        Volume <- pi * (cell.diameter / 2) ^ 2 * cell.initialLength * 2^((t %% mean.generationTime) / mean.generationTime)
        # Drug concentration
        DrugConcentration <- drug / Volume
        
        # Quantities that are repeatedly used
        TargetBinding <- rate.target.binding * target.protein * DrugConcentration
        TargetUnbinding <- rate.dissociation * target.dimer
        TargetDeacylation <- rate.deacylation * target.dimer
        EffluxBinding <- rate.efflux.binding * efflux.protein * DrugConcentration
        EffluxUnbinding <- rate.dissociation * efflux.dimer
        EffluxDeacylation <- rate.deacylation * efflux.dimer
        Catalysis <- rate.catalysis * efflux.dimer
        
        # RHS of differential equations    
        d.target.protein    <- rate.target.protein.translation - TargetBinding + TargetUnbinding + TargetDeacylation
        d.efflux.DNA        <- rate.efflux.DNA.switchOn * (1 - efflux.DNA) - rate.efflux.DNA.switchOff * efflux.DNA
        d.efflux.mRNA       <- rate.efflux.mRNA.transcription * efflux.DNA - rate.efflux.mRNA.deg * efflux.mRNA
        d.efflux.protein    <- rate.efflux.protein.translation * efflux.mRNA - EffluxBinding + EffluxUnbinding + EffluxDeacylation + Catalysis
        d.drug              <- rate.drug.diffusion * SurfaceArea * (ifelse(t < drug.off, drug.outsideConcentration, 0) - DrugConcentration) - (TargetBinding + EffluxBinding) + (TargetUnbinding + EffluxUnbinding)
        d.target.dimer      <- TargetBinding - TargetUnbinding - TargetDeacylation
        d.efflux.dimer      <- EffluxBinding - (EffluxUnbinding + EffluxDeacylation + Catalysis)
        
        return(list(c(d.target.protein, d.efflux.DNA, d.efflux.mRNA, d.efflux.protein, d.drug, d.target.dimer, d.efflux.dimer)))
    })
}

integrateDivisionEqns <- function(t.max, params, IC){
    #' Initial values for each variable
    vars <- c("target.protein", "efflux.DNA", "efflux.mRNA", "efflux.protein", "drug", "target.dimer", "efflux.dimer")
    init.cond <- unlist(IC[1, vars])
    #' Times
    times <- seq.int(0, t.max, by = params$mean.generationTime/100)
    # #' Function for what happens at division time
    # divideFun <- function(t, y, parms){
    #     # Divide all values except DNA by 2
    #     div.vars <- (names(y) != "efflux.DNA")
    #     y[div.vars] <- y[div.vars] / 2
    #     return(y)
    # }
    # # Function to determine division time (excluding division at t=0)
    # rootFun <- function(t, y, parms) {
    #     return( (t %% parms$mean.generationTime) + (t == 0) )
    # }
    # # Solve
    # sol <- as.data.frame(ode(init.cond, times, divisionEqns, params,
    #                          events = list(func = divideFun, root = TRUE),
    #                          rootfun = rootFun,
    #                          hmax = 5e-3, atol = 1e-10)
    #                      )
    
    #' Data frame for division events, which halves the selected variables at each division time
    division.events <- DFapply(seq(from = params$mean.generationTime, to = t.max, by = params$mean.generationTime), function(div.time){
        data.frame(var = vars[vars != "efflux.DNA"],
                   time = div.time,
                   value = 0.5,
                   method = "multiply"
        )
    })
    #' Solve
    sol <- as.data.frame(ode(init.cond, times, divisionEqns, params,
                             events = list(data = division.events))
    )
    #' Add volume
    sol$volume <- with(params, pi * (cell.diameter / 2) ^ 2 * cell.initialLength * 2 ^ ((sol$time %% mean.generationTime) / mean.generationTime))
    
    #' Return
    return(sol)
}


integrateDivisionEqns <- function(t.max, params, IC){
    # Initial values for each variable
    vars <- c("target.protein", "efflux.DNA", "efflux.mRNA", "efflux.protein", "drug", "target.dimer", "efflux.dimer")
    init.cond <- unlist(IC[1, vars])
    # Times
    times <- seq.int(0, t.max, by = params$mean.generationTime/100)

    
    #' Data frame for division events, which halves the selected variables at each division time
    division.events <- DFapply(seq(from = params$mean.generationTime, to = t.max, by = params$mean.generationTime), function(div.time){
        data.frame(var = vars[vars != "efflux.DNA"],
                   time = div.time,
                   value = 0.5,
                   method = "multiply"
                   )
    })
    
    # Solve
    sol <- as.data.frame(ode(init.cond, times, divisionEqns, params,
                             events = list(data = division.events))
    )
    # Add volume
    sol$volume <- with(params, pi * (cell.diameter / 2) ^ 2 * cell.initialLength * 2 ^ (((sol$time - 1e-10) %% mean.generationTime) / mean.generationTime))
    
    # Return
    return(sol)
}
