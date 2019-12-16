#!/bin/bash

homepath=/home/e924p726/My-Code-Collection/

rm -r bin/
mkdir bin/
mkdir bin/wham_gen
mkdir bin/utils
mkdir bin/msd
mkdir bin/build
mkdir bin/widom
mkdir bin/gofr
mkdir bin/gyro
mkdir bin/gendist


# Sets up the symlinks
ln -s $homepath/Util/CP2K/wham_generation/setup_wham.py bin/wham_gen
ln -s $homepath/Util/CP2K/wham_generation/wham_histogram.py bin/wham_gen
ln -s $homepath/Util/CP2K/wham_generation/conv.py bin/wham_gen
ln -s $homepath/Util/block_average/calculate_average.py bin/utils
ln -s $homepath/Util/curve_fit/calculate_fit.py bin/utils
ln -s $homepath/Util/Lammps_Log/read_log.py bin/utils
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/calc_msd.exe bin/msd
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/gen_msd_inps.py bin/msd
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/fit_msds.py bin/msd
ln -s $homepath/Simulation_Analysis/Mean-Squared-Displacement/parse_msds.py bin/msd
ln -s $homepath/Util/general_system/build.py bin/build
ln -s $homepath/Util/general_system/molec_generator.py bin/build
ln -s $homepath/Simulation_Codes/MD_Widom/conv_connectivity.py bin/widom
ln -s $homepath/Simulation_Codes/MD_Widom/gen_ins.py bin/widom
ln -s $homepath/Simulation_Codes/MD_Widom/widom_backbone.py bin/widom
ln -s $homepath/Simulation_Codes/MD_Widom/widom_calculation.py bin/widom
ln -s $homepath/Simulation_Codes/MD_Widom/widom_farmer.py bin/widom
ln -s $homepath/Simulation_Analysis/Pair-Distribution/calc_gofr.exe bin/gofr
ln -s $homepath/Simulation_Analysis/Pair-Distribution/calc_egofr.exe bin/gofr
ln -s $homepath/Simulation_Analysis/Pair-Distribution/calc_e2gofr.exe bin/gofr
ln -s $homepath/Simulation_Analysis/Pair-Distribution/setup_gofr.py bin/gofr
ln -s $homepath/Simulation_Analysis/Pair-Distribution/manipulate_gofr.py bin/gofr
ln -s $homepath/Simulation_Analysis/Atomic-Gyration/calc_atom_gyro.exe bin/gyro
ln -s $homepath/Simulation_Analysis/Atomic-Gyration/setup_gyro.py bin/gyro
ln -s $homepath/Simulation_Analysis/Predict_Distributions/general_distribution.py bin/gendist

chmod 777 bin/msd/*
chmod 777 bin/wham_gen/*
chmod 777 bin/build/*
chmod 777 bin/widom/*
chmod 777 bin/utils/*
chmod 777 bin/gofr/*
chmod 777 bin/gyro/*
chmod 777 bin/gendist/*
