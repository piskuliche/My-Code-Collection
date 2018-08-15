from mpi4py import MPI
from lammps import PyLammps

L = PyLammps()

L.file('in.melt')



L.command('run 10')

if MPI.COMM_WORLD.rank == 0:
    pe = L.eval("pe")
    print("Potential Energy:", pe)

