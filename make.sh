#!/bin/bash

homepath=/home/e924p726/My-Code-Collection/

rm -r bin/
mkdir bin/


# Sets up the symlinks
ln -s $homepath/Util/CP2K/wham_generation/setup_wham.py bin/
ln -s $homepath/Util/CP2K/wham_generation/wham_histogram.py bin/
ln -s $homepath/Util/block_average/calculate_average.py bin/
ln -s $homepath/Util/curve_fit/calculate_fit.py bin/
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/calc_msd.exe bin/
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/gen_msd_inps.py bin/
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/fit_msds.py bin/
ln -s $homepath/Util/general_system/build.py bin/
ln -s $homepath/Simulation_Codes/MD_Widom/conv_connectivity.py bin/
ln -s $homepath/Simulation_Codes/MD_Widom/gen_ins.py bin/
ln -s $homepath/Simulation_Codes/MD_Widom/widom_backbone.py bin/
ln -s $homepath/Simulation_Codes/MD_Widom/widom_calculation.py bin/
ln -s $homepath/Simulation_Codes/MD_Widom/widom_farmer.py bin/
ln -s $homepath/Util/general_system/molec_generator.py bin/

chmod 777 bin/*
