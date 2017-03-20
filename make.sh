#!/bin/bash

#This bash utility is written to simply make sure that the programs within this directory
#end up in the correct folders (based on the path environment variables)

FINITE_PATH=/panfs/pfs.local/work/thompson/e924p726/Programs

echo "Welcome to the Python Code Collection: Unpacking"

# Creates the Python Directory if it doesn't exist.
if [ ! -d "$FINITE_PATH/Python" ]; then
    mkdir $FINITE_PATH/Python
    mkdir $FINITE_PATH/Python/Executable $FINITE_PATH/Python/Not-Executable
    echo "$FINITE_PATH/Python not found - Creating."
fi

# CLEARS THE FOLDER
rm $FINITE_PATH/Python/* $FINITE_PATH/Python/Executable/* $FINITE_PATH/Python/Not-Executable/*

# Copies the programs into the proper directory
# Executable Programs
cp Block-Average/calculate_average.py $FINITE_PATH/Python/Executable
cp MSD-TO-D/MSD_to_D.py $FINITE_PATH/Python/Executable
cp Minimization/minimize.py $FINITE_PATH/Python/Executable

# Not-Executable Programs
cp Curve-Fitting/calculate_fit.py $FINITE_PATH/Python/Not-Executable


# Creates a New Bash Code Chunk
if [ -e "~/.bash_PCC" ]; then
    rm "~/.bash_PCC"
fi

echo 'export PATH="'$FINITE_PATH'/Python/Executable:$PATH"' > ~/.bash_PCC

# Append to BashRC
if ! grep -q PCCFLAG ~/.bashrc ; then
    cat Src/brc_mod >> ~/.bashrc
fi 

echo "Python Code Collection: Unpack Complete"
