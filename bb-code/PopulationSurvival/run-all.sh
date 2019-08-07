#!/bin/bash

echo "WT"
./simulate.sh PopulationSurvival/config-ciprofloxacin-WT.R
echo "REG-BURST"
./simulate.sh PopulationSurvival/config-ciprofloxacin-REG-BURST.R
echo "REG-ON"
./simulate.sh PopulationSurvival/config-ciprofloxacin-REG-ON.R
echo "REG-OFF"
./simulate.sh PopulationSurvival/config-ciprofloxacin-REG-OFF.R
echo "STRUCT-BIND"
./simulate.sh PopulationSurvival/config-ciprofloxacin-STRUCT-BIND.R
echo "STRUCT-CAT"
./simulate.sh PopulationSurvival/config-ciprofloxacin-STRUCT-CAT.R
echo "KO"
./simulate.sh PopulationSurvival/config-ciprofloxacin-KO.R
