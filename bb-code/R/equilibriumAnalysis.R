#' equilibriumAnalysis.R
#' Author: Peter Ashcroft, ETH Zurich

#' Dependencies:


#' Function list:
#' - getEquilibriumValues()


#' Extract the equilibrium values of DNA, mRNA, and Protein for wildtype cells
getEquilibriumValues <- function(model.params){
    # TARGET
    #target.DNA.eq <- with(model.params, rate.target.DNA.switchOn / (rate.target.DNA.switchOn + rate.target.DNA.switchOff) )
    #target.mRNA.eq <- with(model.params, (rate.target.mRNA.transcription / rate.target.mRNA.deg) * target.DNA.eq)
    #target.protein.eq <- with(model.params, (rate.target.protein.translation / rate.target.protein.deg) * target.mRNA.eq)
    target.protein.eq <- model.params$target.protein.total
    # EFFLUX
    efflux.DNA.eq <- with(model.params, rate.efflux.DNA.switchOn / (rate.efflux.DNA.switchOn + rate.efflux.DNA.switchOff) )
    efflux.mRNA.eq <- with(model.params, (rate.efflux.mRNA.transcription / (rate.efflux.mRNA.deg + log(2)/generationTime) ) * efflux.DNA.eq)
    efflux.protein.eq <- with(model.params, (rate.efflux.protein.translation / (rate.efflux.protein.deg + log(2)/generationTime) ) * efflux.mRNA.eq)
    
    # DRUG
    concentration.eq <- with(model.params, drug.outsideConcentration / (1 + (1/rate.drug.diffusion) * (log(2)/generationTime) * (cell.diameter/4) ) )
    drug.eq <- with(model.params, concentration.eq * pi * ((cell.diameter^2)/4) * (2/log(2))) # 2/log(2) is the mean cell length in \mu m
    
    # TARGET--DRUG dimer
    #target.dimer.eq <- with(model.params, (rate.target.dimer.creation/rate.target.dimer.dissociation) * target.protein.eq * concentration.eq )
    
    # Combine to dataframe
    mol.names <- c("target.protein", "efflux.DNA", "efflux.mRNA", "efflux.protein", "drug", "concentration")
    mol.number <- c(target.protein.eq, efflux.DNA.eq, efflux.mRNA.eq, efflux.protein.eq, drug.eq, concentration.eq)
    equilibrium.values <- data.frame(molecule = mol.names, number = mol.number, stringsAsFactors = FALSE)
    
    return(equilibrium.values)
}
