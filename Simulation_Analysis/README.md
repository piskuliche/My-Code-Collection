# Simulation Analysis Codes:

In this directory there are multiple codes that can be used with a molecular simulation to calculate important quantities.

## Atomic-Gyration

This code reads a LAMMPS trajectory and calculates the simulation box radius of gyration.

## Mean-Squared-Displacement

This code calculates the MSD of an arbitrary set of molecules from a lammps trajectory file, and calculates the diffusion coefficient of those molecules.

## Pair-Distribution

This code calculates the radial distribution function (and its temperature,pressure derivatives) from a lammps simulation using a trajectory file and the lammps log file.

## Predict-Distributions

This code calculates the Roo, Roh, Ctheta, and Asphericity distributions for water Hbond angles, as well as the radial distribution function. It also calculates the first derivative of all of these quantities.

## VMD-SCRIPTS

This is a collection of VMD tcl scripts that can be used with simulations.

## Viscosity

This code uses the lammps log file to calculate the shear viscosity using a green-kubo relation

## Minimization

This code calculates the minimum along a particular axis of a vector.

