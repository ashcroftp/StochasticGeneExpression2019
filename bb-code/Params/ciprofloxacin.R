#' ciprofloxacin.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' Parameters that define the action of the drug.
#' Unit conversions are important! Units are given in square brackets,
#' and are standardised into number of molecules [n], grams [g], micrometers [um], minutes [min].
#' 
#' Values from the literature ----
#' Number of targets:       ~300 DNA-gyrase complexes per cell [n] [1, 2]
#' MIC drug concentration:  0.03 [mg L^{-1}] [3]
#' Drug molar mass:         331.347 [g mol^{-1}] [4]
#' Diffusion rate:          2.0*10^{-11} [m s^{-1}] [5]
#' Binding rate:            ~3.6*10^{4} [M^{-1} s^{-1}] [6: Based on the cited paper's derivations, binding constant Kq is similar to the inhibition constant Kiq = 8.4 [nM]; binding rate can be obtained by dividing the dissociation rate (shown in next line) with Kiq]
#' Dissociation rate:       3*10^-4 [s^{-1}] [6: Experimentally determined]
#' Catalysis rate:          10 [s^{-1}] [7: This value is based on AcrAB efflux pump that confers multi-drug resistance]
#' Deacylation rate:        0 [s^{-1}] [N/A]
#' Minimum growth rate:     -6.5 [hr^{-1}] [3: When under very high drug concentration, measured as the change of bacterial density per hour in log 10]
#' 
#' Conversion formulae ----
#' Targets:                 We use N_0 = (2/3)*N such that N(0)=N_0, N(t_G)=2N_0, and N(t_G/2)=N.
#' Concentration:           1 [mg L^{-1}] = 10^{-18} [g um^{-3}]
#' Molar mass:              1 [g mol^{-1}] = 1/N_A [g n^{-1}]
#' Molecular concentration: concentration / molar mass [n um^{-3}] (10^{-18}*N_A = 6.022*10^5)
#' Diffusion rate:          1 [m s^{-1}] = 60 * 10^6 [um min^{-1}]
#' Binding rate:            1 [M^{-1} s^{-1}] = 60 * 10^15 / N_A [um^3 n^{-1} min^{-1}]
#' Dissociation rate:       1 [s^{-1}] = 60 [min^{-1}]
#' Growth rate:             1 [hr^{-1}] = 1/60 [min^{-1}]
#' 
#' References ----
#' [1] -- Chong et al., Cell, 158(2), 314-326 (2014)
#' [2] -- Stracy et al., bioRxiv (https://doi.org/10.1101/460006), (2018)
#' [3] -- Regoes et al., Antimicrob. Agents Chemother., 48(10), 3670-3676 (2004)
#' [4] -- National Center for Biotechnology Information. PubChem Compound Database; CID=2764, https://pubchem.ncbi.nlm.nih.gov/compound/2764 (accessed May 16, 2018).
#' [5] -- Abel zur Wiesch et al., Sci. Transl. Med. 7, 287ra73 (2015)
#' [6] -- Kampranis & Maxwell, J. Biol. Chem., 273(35), 22615-22626 (1998)
#' [7] -- Nagano & Nikaido, PNAS, 106(14), 5854-5858 (2009)
#' 

#' Computed values to be loaded into namespace
#' Number of targets
drug.targets <- round(300 * (2/3)) # [n]
#' Drug effect
drug.MIC <- round((0.03 * 1.0e-18) / (331.347 / 6.022e23)) # [n um^{-3}]: MIC drug concentration
drug.kappa <- 3.0 # [dimensionless]: Shape (steepness) parameter (\kappa)
#' Drug dynamics
drug.outsideConcentration <- 0.0 # [n um^{-3}]: Concentration of drug in the environment (constant: c_out, zero by default)
drug.diffusionConstant <- 2.0e-11 * 60.0*1.0e6 # [um minute^{-1}]: Rate at which drug enters the cell (\sigma)
#' Dimer parameters
drug.bindingRate <- 3.6e4 * 60.0*1.0e15/6.022e23  # [um^3 n^{-1} minute^{-1}]: Rate protein binds to drug (k_f)
dimer.dissociationRate <- 3.0e-4 * 60.0 # [minute^{-1}]: Dissociation of dimer into protein and drug (k_b)
dimer.catalysisRate <- 10 * 60.0 # [minute^{-1}]: Removal/catalysis of the drug by an enzyme (k_cat)
drug.deacylationRate <- 0 # [minute^{-1}]: Drug will spontaneously de-activate and release an unbound target or efflux (k_d)
#' Population dynamics
pop.minGrowthRate <- -6.5 / 60.0 # [minute^{-1}]: Minimum growth rate (\psi_{min})
#' MIC fraction of bound targets
MIC.fraction <- 0.085 # Determined by simulation (\rho_{MIC})

#efflux.binding.factor <- 2
