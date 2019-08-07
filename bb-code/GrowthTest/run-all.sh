#!/bin/bash

echo "WT"
./simulate.sh GrowthTest/config-WT-growth.R
echo "REG-BURST"
./simulate.sh GrowthTest/config-REG-BURST-growth.R
echo "REG-ON"
./simulate.sh GrowthTest/config-REG-ON-growth.R
echo "REG-OFF"
./simulate.sh GrowthTest/config-REG-OFF-growth.R
echo "STRUCT-BIND"
./simulate.sh GrowthTest/config-STRUCT-BIND-growth.R
echo "STRUCT-CAT"
./simulate.sh GrowthTest/config-STRUCT-CAT-growth.R
echo "KO"
./simulate.sh GrowthTest/config-KO-growth.R
