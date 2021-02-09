### Occupancy-Count

This is a program that is intended to calculate the occupancy of water molecules surrounding some sort of solute molecule.

It was initially intended to work with Acrylimide polymers, but could be relatively simply be modified to other uses.

It takes two input files: 

1) inp_vals
2) solv.in

As well as a trajectory file.

## inp_vals

This file is formatted as follows:
1st line: number of atoms in solute
2nd-nth lines: one entry per line of atom name and distance defining cutoff i.e. C1 4.5

We have included an example in the examples folder. 

## solv.in

We have included an example input file in the examples folder. 

## Other Code

We have also include a python code that takes the data file from the f90 program and outputs it in a useful manner.


