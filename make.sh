#!/bin/bash

homepath=/home/e924p726/My-Code-Collection/

rm -ri bin/
mkdir bin/


# Sets up the symlinks
ln -s $homepath/CP2K_Util/wham_generation/setup_wham.py bin/
ln -s $homepath/CP2K_Util/wham_generation/wham_histogram.py bin/
ln -s $homepath/block_average/calculate_average.py bin/
ln -s $homepath/curve_fit/calculate_fit.py bin/
ln -s $homepath/system_builds/water/build_water.py bin/
ln -s $homepath/system_builds/acn/build_acn.py bin/
ln -s $homepath/system_builds/acn/cxl/build_acn_cxl.py bin/
ln -s $homepath/system_builds/zeolites/zeolite_xyz_to_lmps.py bin/c

chmod 777 bin/*
