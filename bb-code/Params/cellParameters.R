#' cellParameters.R
#' Authors: Lei Sun & Peter Ashcroft, ETH Zurich

#' Parameters that define the dynamics of the cell.
#' Unit conversions are important! Units are given in square brackets,
#' and are standardised into number of molecules [n], grams [g], micrometers [um], minutes [min].
#' 
#' Values from the literature ----
#' Population growth rate:  0.80 [hr^{-1} (log-10)] [1] # weighted.mean(x = c(0.88, 0.75, 0.70, 0.89, 0.81), w = 1/c(0.16, 0.21, 0.23, 0.41, 0.42))
#' Generation time:         0.37 [hr] [1] # log(2)/(log(10)*weighted.mean(x = c(0.88, 0.75, 0.70, 0.89, 0.81), w = 1/c(0.16, 0.21, 0.23, 0.41, 0.42)))
#' Diameter:                934±6 [nm] [2; Fig. 2A]
#' Initial length:          2.0-4.0 [um] [3; Figs. 2D,E & 5A] # Media-dependent: 2.5-4.0 um in LB, 2.0-3.5 um in M9
#' DNA ON-state:            5.0-60 [min] [4; Table 1] [9; Fig. 2]
#' DNA OFF-state:           0.5-3000 [min] [4; Table 1] [9; Fig. 2]
#' mRNA transcription rate: 1/2.5 [min^{-1}] [5; Fig. 3]
#' mRNA half-life:          3-8 [min] [6; Fig. 1]
#' or, mRNA half-life:      2 [min] [7; Fig. 1]
#' Protein burst size:      5-40 [n] [7; Fig. 1] # Molecules per transcriptional burst
#' Protein degradation:     Protein degradation is negligible compared to dilution-by-cell-division [8]
#' 
#' References ----
#' [1] -- Regoes et al., Antimicrob. Agents Chemother., 48(10), 3670-3676 (2004)
#' [2] -- Ouzounov et al., Biophys. J., 111(5), 1035-43 (2016)
#' [3] -- Campos et al., Cell, 159, 1433-1446 (2004)
#' [4] -- Lionnet & Singer, EMBO reports, 13(4), 313-321 (2012)
#' [5] -- Golding et al., Cell, 123, 1025-1036 (2005)
#' [6] -- Bernstein et al., PNAS, 99(15), 9697-9702 (2002)
#' [7] -- Thattai & van Oudenaarden, PNAS, 98(15), 8614-8619 (2001)
#' [8] -- Larrabee et al., J. Biol. Chem., 255(9), 4125-4130 (1980)
#' [9] -- Hammar et al., Nat. Genet. 46(4), 405-408 (2014)
#' 

#' Generation time, cell size, and growth rate parameters
pop.maxGrowthRate <- weighted.mean(x = c(0.88, 0.75, 0.70, 0.89, 0.81), w = 1/c(0.16, 0.21, 0.23, 0.41, 0.42)) / 60.0 # [min^{-1}] (\psi_{min})
cell.generationTime <- round(log(2)/(log(10)*pop.maxGrowthRate), digits = 1) # [min]: Mean time between cell divisions (rounded for simmulation reasons) (t_G)
cell.diameter <- 934 / 1e3 # [um]: diameter of cylindrical cell (2*r)
cell.initialLength <- 3.4 # [um]: length of cell after division [3; Fig. 2E] (\ell_0)
#' Parameters of EFFLUX DNA, mRNA, and protein
#' DNA parameters
#' Transcription lasts, on average, 5 min, with 195 min OFF-interval. 
#' This means transcription is ON 5/(5+195)=2.5% of all time.
efflux.DNA.onTime <- 5.0 # [min]
efflux.DNA.offTime <- 195.0 # [min]
efflux.DNA.switchOff <- 1.0 / efflux.DNA.onTime # [min^{-1}] (\Gamma_{off})
efflux.DNA.switchOn <- 1.0 / efflux.DNA.offTime # [min^{-1}] (\Gamma_{on})
#' mRNA parameters
efflux.mRNA.transcriptionRate <- (1/2.5) # [min^{-1}] (\tau)
efflux.mRNA.halfLife <- 5.0 # [min] From [6] 
#efflux.mRNA.halfLife <- 2.0 # [min] From [7]
efflux.mRNA.degRate <- log(2) / efflux.mRNA.halfLife # [min^{-1}] (\gamma)
#' Protein parameters
efflux.protein.burstSize <- (5 + 40)/2 # [n] (b)
efflux.protein.translationRate <- efflux.protein.burstSize * efflux.mRNA.degRate # [n min^{-1}] Derivation from [7] (\gamma b)
#efflux.protein.halfLife <- 60.0 # [min] From [7]
#efflux.protein.degRate <- log(2) / efflux.protein.halfLife # [min^{-1}]
efflux.protein.degRate <- 0 # [min^{-1}] Neglecting degradation as it is insignificant compared to dilution-by-cell-division, based on [8]
#' Multiplicative effect of mutants
mutant.effect <- 2.0 # [dimensionless] (\mu)
#' Fitness costs of mutants (percentage)
fitness.cost <- 0.05 # [dimensionless] (\nu)

#' Still not sure that this is needed...
efflux.binding.factor <- 1
