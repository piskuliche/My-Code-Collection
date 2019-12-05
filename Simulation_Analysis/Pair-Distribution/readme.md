# README

## Compilation

To compile this code, use the following:

```
pgfortran -O3 -mcmodel=medium epairdist.f90 -o calc_gofr.exe
```


## Instructions

To run this code you need four files (for a gofr calculation) or five files (for a derivative calculation).

The files always required are:
gofr.inp: simulation input file in a namelist format, can be created by setup_gofr.py
traj.file: trajectory file in an xyz format
L.dat: File that specifies the box length at each step (1 length per frame per line)
molinfo.dat: special file created by setup_gofr.py from the data file, basically just the atom information is included.

The files that are only required in the case of a derivative calculation is:
X_init.out where X is the type of energy that is being used for the derivative.


