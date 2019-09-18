#!/usr/bin/env bash

# Set library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ajasper/lib/

# Set onedmin.x executable path
ONEDMINEXE=/lcrc/project/CMRP/pacc/OneDMin/build/onedmin.x

# Run several onedmin.x instances
cd /lcrc/project/CMRP/pacc/OneDMin/examples/CH4_N2/parallel/run/SPC/CH4/VNWKTOKETHGBQD/0/1/UHFFFAOYSA-N/ETRANS/3vx11mR/run1
time $ONEDMINEXE < input.dat > output.dat &
cd ../run2
time $ONEDMINEXE < input.dat > output.dat &
cd ../run3
time $ONEDMINEXE < input.dat > output.dat &
cd ../run4
time $ONEDMINEXE < input.dat > output.dat &
wait

