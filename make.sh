#!/bin/bash

homepath=/home/e924p726/My-Code-Collection/

rm -r bin/
mkdir bin/


# Sets up the symlinks
ln -s $homepath/CP2K_Util/wham_generation/setup_wham.py bin/
ln -s $homepath/CP2K_Util/wham_generation/wham_histogram.py bin/
ln -s $homepath/block_average/calculate_average.py bin/
ln -s $homepath/curve_fit/calculate_fit.py bin/
ln -s $homepath/system_builds/water/build_water.py bin/
ln -s $homepath/system_builds/acn/build_acn.py bin/
ln -s $homepath/system_builds/acn/cxl/build_acn_cxl.py bin/
ln -s $homepath/system_builds/zeolites/zeolite_xyz_to_lmps.py bin/
ln -s $homepath/system_builds/acn/cxl/withions/build_acn_cxl_ion.py bin/
ln -s $homepath/diffusion/msd/calc_msd.exe bin/
ln -s $homepath/diffusion/msd/gen_msd_inps.py bin/
ln -s $homepath/diffusion/msd/fit_msds.py bin/
ln -s $homepath/system_builds/general_system/build.py bin/
ln -s $homepath/MD_Widom/conv_connectivity.py bin/
ln -s $homepath/MD_Widom/gen_ins.py bin/
ln -s $homepath/MD_Widom/widom_backbone.py bin/
ln -s $homepath/MD_Widom/widom_calculation.py bin/
ln -s $homepath/system_builds/general_system/molec_generator.py bin/

chmod 777 bin/*
