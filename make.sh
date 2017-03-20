#!/bin/bash

#This bash utility is written to simply make sure that the programs within this directory
#end up in the correct folders (based on the path environment variables)

FINITE_PATH=/panfs/pfs.local/work/thompson/e924p726/Programs

# Creates the Python Directory if it doesn't exist.
if [ ! -d "$FINITE_PATH/Python" ]; then
    mkdir $FINITE_PATH/Python
    mkdir $FINITE_PATH/Python/Executable $FINITE_PATH/Python/Not-Executable
fi

# CLEARS THE FOLDER
rm $FINITE_PATH/Python/* $FINITE_PATH/Python/Executable/* $FINITE_PATH/Python/Not-Executable/*

# Copies the programs into the proper directory
# Executable Programs
cp Block-Average/calculate_average.py $FINITE_PATH/Python/Executable
cp 
# Not-Executable Programs
cp Curve-Fitting/calculate_fit.py $FINITE_PATH/Python/Not-Executable
