#!/bin/bash

echo "WT"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-WT-growth.R
echo "REG-BURST"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-REG-BURST-growth.R
echo "REG-ON"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-REG-ON-growth.R
echo "REG-OFF"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-REG-OFF-growth.R
echo "STRUCT-BIND"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-STRUCT-BIND-growth.R
echo "STRUCT-CAT"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-STRUCT-CAT-growth.R
echo "KO"
R --vanilla --slave < ../Simulate/simulate.R --args PopulationGrowth/config-KO-growth.R
