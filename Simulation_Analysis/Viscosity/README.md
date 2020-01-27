# Shear Viscosity Correlation Function

```python
python viscosity.py -f log.lammps -nb 5 ...
```

This code takes a lammps log file and extracts the data from it for the pressure tensor and
uses it to calculate the shear viscosity Green-Kubo relation.

This produces two output file sets, one is shear.dat which stores the average and uncertainty as well as the unintegrated form.

The second is a blocked version of the shear function.
