#!/bin/bash

echo "WT"
./simulate.sh PopulationGrowth/config-ciprofloxacin-WT.R
echo "REG-BURST"
./simulate.sh PopulationGrowth/config-ciprofloxacin-REG-BURST.R
echo "REG-ON"
./simulate.sh PopulationGrowth/config-ciprofloxacin-REG-ON.R
echo "REG-OFF"
./simulate.sh PopulationGrowth/config-ciprofloxacin-REG-OFF.R
echo "STRUCT-BIND"
./simulate.sh PopulationGrowth/config-ciprofloxacin-STRUCT-BIND.R
echo "STRUCT-CAT"
./simulate.sh PopulationGrowth/config-ciprofloxacin-STRUCT-CAT.R
echo "KO"
./simulate.sh PopulationGrowth/config-ciprofloxacin-KO.R
