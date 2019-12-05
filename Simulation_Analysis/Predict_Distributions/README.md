### Distribution Prediction Software

## Introduction

In brief, this code calculates a probability distribution (P) and calculates its first derivative with
respect to temperature (or technically pressure).

Currently, this code supports (and calculates) the following distributions:
1) OO Hydrogen Bond Donor-Acceptor Distance Distributions
2) HO Hydrogen Bond Donor-Acceptor Distance Distributions
3) OHO Hydrogen Bond Angle Distributions
4) Voronoi Asphericity (adds an order of magnitude to calculation time)

The code, in general, is written to be modular. For properties that primarily depend on positional data it should be easy to add further calculations and achieve success.

## Code Dependencies

There are a couple of files that the code expects to do its job, these are:
1) e\_init.out (file with 1 energy per frame per line)
2) trajectory file
3) L.dat (file with box lengths in the same format as e\_init.out)

There are also a few optional files:
1) lj\_init.out (file with Lennard-Jones energies)
2) ke\_init.out (file with Kinetic energies)
3) vol\_init.out (file with Volumes)

## Arguments 
Arguments will be listed when the code is run with the optional -h argument.

All arguments include defaults; however, for best results one should carefully select arguments.

**-h**: prints help information

**-f**: trajectory file name

**-nblocks**: number of blocks for block averaging

**-nconfigs**: number of configurations

**-oatm**: integer type representing oxygen in the trajectory file

**-hatm**: integer type representing hydrogen in the trajectory file

**-order**: order of the trajectory file (currently, only 'ohh' is supported and is default)

**-T**: Temperature of simulation, used in calculation of free energy

**-P**: Pressure of simulation, currently unused.

**-prepend**: String to prepend to output files.

**-restart**: 0 if not restarting a calcualtion, 1 if restarting a farmed calculation, 2 if launching a farmed calculation

**-restno**: integer value that is the number of farmed trajectories.

**-rest\_freq**: integer value, the number of configurations between dumping restart files (these overwrite each other)

**-frame\_freq**: integer value, the number of configurations per subdirectory


## Running the Code (Normal Run):

For calculations that are simple, i.e. calculations that don't take long to converge, it is very simple to call the code. 



For a simple run that includes asphericity (and default hbond distributions)
general\_distribution.py -f "traj.xyz" -nblocks 5 -nconfigs 1000 -prepend run\_ -asphere 1

The code will run, and all output files will be prepended with run\_.


## Running the Code (Farming Run):

For calculations that are complicated or take a long time to run, the code has a built in function to speed up the 
calculation that relies on running the code three separate times.

1) Farm Launch
The first run of the code aims to split the trajectory file into separate files by creating a subdirectory Farm with frame_freq configs per subdirectory.

i.e. this command,
general\_distribution.py -f traj.xyz -nblocks 5 -nconfigs 1000000 -prepend farm1\_ -restart 2 -frame\_freq 1000 

will separate the trajectory file, box length file, and energy file into 1000 directories each that has 1000 configurations.


2) Calculate distributions in subdirectories
This is best accomplished by a job array, where you individually run the code in each subdirectory.

i.e. for the above example

cd Farm/$SLURM\_ARRAY\_TASK\_ID

general\_distribution.py -f traj.xyz -nblocks 5 -nconfigs 1000 -prepend farm1\_ -asphere 1 -rest\_freq 1000 > out

cd ../../

Note that nconfigs has changed to the number of configurations in the particular subdirectory.

3) Read the Restarts and Combine
The last run of the code reconstitutes the separate calculations and reads in their results and finalizes the calculations.

general\_distribution.py -f traj.xyz -nblocks 5 -T 298.15 -nconfigs 1000000 -prepend farm1\_ -asphere 1 -restart 1 -restno 1000

Note that the nconfigs has returned to the total amount. This code should be run *inside* the farm directory.



