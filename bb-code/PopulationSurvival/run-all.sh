#!/bin/bash

R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-WT.R
R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-KO.R
R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-REG-ON.R
R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-REG-OFF.R
R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-REG-BURST.R
R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-STRUCT-CAT.R
R --vanilla --slave < ../Simulate/simulate.R --args PopulationSurvival/config-ciprofloxacin-STRUCT-BIND.R

