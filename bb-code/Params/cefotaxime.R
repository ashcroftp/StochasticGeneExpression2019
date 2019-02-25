#' cefotaxime.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' Parameters that define the action of the drug.
#' Unit conversions are important! Units are given in square brackets,
#' and are standardised into number of molecules [n], grams [g], micrometers [um], minutes [min].
#' 
#' Values from the literature ----
#' Number of targets:       ~200 PBP3 per cell [n] 
#' [1: 200 is the reported value in this paper; it also mentioned two earlier studies with slightly lower values; we use 200 since it's the most recent, and its reported cell growth condition is closest to our simulations]
#' [2: CTX has highest selectivity for PBP3, as shown in Table 2]
#' MIC drug concentration:  ~0.06 [mg L^{-1}] 
#' [3: a wide range of MIC values have been reported from 0.015 (Ref 2) to over 0.5 (Ref 3); the median MIC however is around 0.06, as shown in Fig 1c of Ref 3]
#' [4: Fig 2a; a more recent paper corroborating Ref 3]
#' Drug molar mass:         455.46 [g mol^{-1}] [5]
#' Diffusion rate:          6.0*10^{-8} [m s^{-1}] [6: an estimated general diffusion rate through the outermembrane of E.coli cell]
#' Binding rate:            1.15-8 *10^{4} [M^{-1} s^{-1}] [7: two values are reported, 8e4 is based on a truncated PBP3, 1.15e4 based on wt PBP3 but at 30 degree Celsius; we take their mean value, which is about 4.6]
#' Dissociation rate:       0 [s^{-1}] [11: Other PBP-targetting drugs bind irreversibly, so we assume the same is true for cefotaxime]
#' Catalysis rate:          2.5 [s^{-1}] [8: This value is based on the TEM-1 beta-lactamase enzyme]
#' Deacylation rate:        2.5*10^{-4} [s^{-1}] [9: Table 1; this value is based on E.coli and PBP1b, the target with 2nd highest specificity for CTX after PBP3 (see Ref 2); the deacylation rate of CTX-PBP3 has not been determined, but should be similarly low]
#' Minimum growth rate:     approx. -4.0 [hr^{-1}] [10: Figure 1 compares well with Figure 2B of Regoes et al. Antimicrob. Agents. Chemother. (2004), where CFU decreases from 10^6 to 10^4 in the first hour. Therefore we use the value extrapolated by Regoes et al.]
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
#' [1] -- Weiss et al., Mol. Microbiol., 25(4), 671-681 (1997)
#' [2] -- Kocaoglu & Carlson, Antimicrob. Agents Chemother., 59(5), 2785-2790 (2015)
#' [3] -- Nath et al., Can. J. Infect. Dis., 6(1) 21-27 (1995)
#' [4] -- Kronvall, J. Clin. Microbiol., 48(12), 4445-4452 (2010)
#' [5] -- National Center for Biotechnology Information. PubChem Compound Database; CID=5742673, https://pubchem.ncbi.nlm.nih.gov/compound/5742673 (accessed June 5, 2018).
#' [6] -- Nagano & Nikaido, PNAS, 106(14), 5854-5858 (2009)
#' [7] -- Adam et al., J. Bacteriol., 179(19), 6005-6009 (1997)
#' [8] -- Saves et al., Biochemistry, 34(37), 11660-11667 (1995)
#' [9] -- Terrak et al., Mol. Microbiol., 34(2), 350-364 (1999)
#' [10] -- Xuan et al., Antimicrob. Agents Chemother., 41(7), 1512-1516 (1997)
#' [11] -- Abel zur Wiesch et al., Sci. Transl. Med. 7, 287ra73 (2015)

#' Computed values to be loaded into namespace
#' Number of targets
drug.targets <- round(200 * (2/3)) # [n]
#' Drug effect
drug.MIC <- round((0.06 * 1.0e-18) / (455.46 / 6.022e23)) # [n um^{-3}]: MIC drug concentration
drug.kappa <- 0.0 # [dimensionless]: Shape (steepness) parameter (\kappa)
#' Drug dynamics
drug.outsideConcentration <- 0.0 # [n um^{-3}]: Concentration of drug in the environment (constant: c_out, zero by default)
drug.diffusionConstant <- 6.0e-8 * 60.0*1.0e6 # [um minute^{-1}]: Rate at which drug enters the cell (\sigma)
#' Dimer parameters
drug.bindingRate <- 4.6e4 * 60.0*1.0e15/6.022e23  # [um^3 n^{-1} minute^{-1}]: Rate protein binds to drug (k_f)
dimer.dissociationRate <- 0 * 60.0 # [minute^{-1}]: Dissociation of dimer into protein and drug (k_b)
dimer.catalysisRate <- 2.5 * 60.0 # [minute^{-1}]: Removal/catalysis of the drug by an enzyme (k_cat)
drug.deacylationRate <- 2.5e-4 * 60.0 # [minute^{-1}]: Drug will spontaneously de-activate and release an unbound target or efflux (k_d)
#' Population dynamics
pop.minGrowthRate <- -4.0 / 60.0 # [minute^{-1}]: Minimum growth rate (\psi_{min})
#' MIC fraction of bound targets
MIC.fraction <- 0.93 # Determined by simulation (\rho_{MIC})

#efflux.binding.factor <- 2