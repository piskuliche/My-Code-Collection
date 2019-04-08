README - CP2K

This is hopefully going to turn into a guide to running CP2K MD simulations.


---------------------------------
---------------------------------
Basic Sections MD:
---------------------------------
---------------------------------

The FORCE_EVAL section is where most of the force calculation parameters are. An example is included below for SPCE water. This will be a commented input file.

&FORCE_EVAL
  METHOD FIST
  &MM 
    &FORCEFIELD 						! Defines the forcefield parameters.
     &SPLINE				
       EMAX_SPLINE 500					! Specify the maximum value of the potential up to which splines will be constructed 
     &END SPLINE
     &BEND								! Specify bonds even if constrained. If constrained - set force constant to zero
        ATOMS H O H
        K 0.
        THETA0 1.8
      &END BEND
      &BOND
        ATOMS O H
        K 0.
        R0 1.8
      &END BOND
      &CHARGE							! These are also defined in the non-bonded file - so not a huge deal.
        ATOM O
        CHARGE -0.8476
      &END CHARGE
      &CHARGE
        ATOM H
        CHARGE 0.4238
      &END CHARGE
      &NONBONDED						! Can generate these with an outside script.
        @INCLUDE ./genpot_spce
      &END NONBONDED
    &END FORCEFIELD
    &POISSON							! Section for Ewald Parameters
      &EWALD
        EWALD_TYPE SPME					! Single Particle Mesh Ewald
        ALPHA .3						! Alpha*rcut = 3.5, Alpha*L = 7 for 14 Ang alpha of 0.3 is about right
        GMAX 30							! Number of grid points to calculate ewald sum.
        O_SPLINE 6						! order of the beta-Euler spline
      &END EWALD
    &END POISSON
  &END MM 
  &SUBSYS											! Section for coordinates, constraints etc.
    &TOPOLOGY										! Topology lets you read in each type of file
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ./spce_equil.xyz
      CONN_FILE_FORMAT MOL_SET						
      &MOL_SET
        &MOLECULE
         NMOL 96
         CONN_FILE_FORMAT PSF						! PSF Files are the typical way of providing connectivity.
         CONN_FILE_NAME ./topology_fist_WAT.psf
        &END
      &END
    &END TOPOLOGY
    &CELL
      ABC 14.219066 14.219066 14.219066				! Cell dimensions.
    &END CELL
     &KIND H                           				! The KIND sections are useful for AIMD calculations - They are not needed for a classical simulation
       BASIS_SET DZVP-GTH                 			! 		and indeed are not used otherwise.
       POTENTIAL GTH-PBE-q1                 	
     &END KIND                        
     &KIND O                           
       BASIS_SET DZVP-GTH                 
       POTENTIAL GTH-PBE-q6                 
     &END KIND                         
   &END SUBSYS                        
    &PRINT											! Prints to the log file every 100 steps.
      &PROGRAM_RUN_INFO
        &EACH
          MD 100
        &END EACH
      &END PROGRAM_RUN_INFO
    &END PRINT
 &END FORCE_EVAL


The GLOBAL section defines base run parameters.

&GLOBAL
  PROJECT free										! Project name
  WALLTIME 23:50:00									! Walltime 
  RUN_TYPE MD  										! Type of run 
  PRINT_LEVEL LOW									! Information printed out 
&END GLOBAL

The MOTION section defines the MD ensemble and other time integration fixes.

&MOTION
  &CONSTRAINT										! Defines the shake constraint
    CONSTRAINT_INIT T								! Turns the constraint on or off.
    SHAKE_TOLERANCE 1E-10							! Tolerance of shake - 1E-10 is probably good.
    &G3X3											! This section defines the collective variables that hold OH, HH bonds rigid.
      DISTANCES 1.8897268 1.8897268 3.0859239
      MOLECULE 1
      ATOMS 1 2 3
    &END G3X3
  &END CONSTRAINT
  &MD
    ENSEMBLE NVT									! Defines the ensemble NVE, NPT, NVT
    STEPS 50000000								 	! Number of steps.
    TIMESTEP 0.5									! Timestep of simulation
    TEMPERATURE 300									! Temperature of the simulation
    !COMVEL_TOL 1.E-7
    &THERMOSTAT										! Defines the thermostat parameters.
      REGION GLOBAL									! The region that the thermostat applies to (whether it is on every molecule or on the whole sys)
      TYPE NOSE										! Defines the Nose hoover thermostat
      &NOSE
        LENGTH 3									! Length of the nose hoover chain.
        YOSHIDA 3									! The order of the yoshida integrator used.
        MTS 2										! Number of multiple timesteps used for the thermostat.
        TIMECON [wavenumber_t] 1000 				! Time constant of the thermostat. 
      &END NOSE
    &END THERMOSTAT
    &PRINT											! Prints integration info every 1000, and energy every 1000.
      &PROGRAM_RUN_INFO
        &EACH
          MD 1000
        &END EACH
      &END PROGRAM_RUN_INFO
      &ENERGY
        &EACH
          MD 1000
        &END EACH
      &END ENERGY
    &END PRINT
  &END MD
  &PRINT											! Print statements for Velocities and Trajectory.
    &VELOCITIES
      &EACH
        MD 500
      &END EACH
    &END VELOCITIES
    &TRAJECTORY
      &EACH
        MD 500
      &END EACH
    &END TRAJECTORY
    &RESTART_HISTORY OFF							! If on - saves restart files
    &END RESTART_HISTORY
    &RESTART ON										! Every 10000 steps saves a restart file.
      &EACH
        MD 10000
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION

The last section EXT_RESTART is only used if you are going to restart the simulation.

&EXT_RESTART
  RESTART_FILE_NAME ./restart
&END EXT_RESTART

---------------------------------
---------------------------------
Generate Non-bonded Interaction File
---------------------------------
---------------------------------

Below is a program that can create the genpot file needed to run cp2k. You can run it like python genpot_cp2k.py

####################################################################################################
# This is a code to create the genpot section of the CP2K input file
# It needs a simple file called param.in that includes the number of species and parameters.
#i.e.  the file should look like:
#---
#   4    nspecies
#   O   -0.8476   0.1553     3.166
#   H    0.4238   0.0000     0.000
#   Li   1.0000   0.3367     1.582
#   F   -1.0000   0.0074     4.514
#---
#If you do this - then you get a genpot section that you can pipe to an output file of your desire.
# Copyright Zeke Piskulich, University of Kansas, July 2018.
# You can contact me at piskuliche@ku.edu
####################################################################################################
import numpy as np

class param_init:
    def __init__(self, nsp):
        self.type = np.chararray((nsp,nsp),itemsize=2)
        self.type[:] = '  '
        self.q = np.zeros((nsp,nsp))
        self.eps = np.zeros((nsp,nsp))
        self.rm = np.zeros((nsp,nsp))

    def add_param(self, type, q, eps, rm, i):
        self.type[i,i] = type
        self.q[i,i] = q
        self.eps[i,i] = eps
        self.rm[i,i] = rm

    def mix_param(self,nsp):
        for i in range(nsp):
            for j in range(nsp):
                self.eps[i,j]=np.sqrt(self.eps[i,i]*self.eps[j,j])
                self.rm[i,j]=0.5*(self.rm[i,i]+self.rm[j,j])
    
    def print_genpot(self, nsp, rcut):
        for i in range(nsp):
            j=nsp-1
            while j >= i:
                print("\t&GENPOT")
                print("\t ATOMS %s %s" % (self.type[i,i],self.type[j,j]))
                print("\t FUNCTION D0*((R0/R)^12-2.0*(R0/R)^6")
                print("\t PARAMETERS D0 R0")
                print("\t UNITS kcalmol angstrom")
                print("\t VARIABLES R")
                print("\t VALUES %s %s" % (self.eps[i,j], self.rm[i,j]))
                print("\t RCUT %s" % rcut)
                print("\t&END GENPOT")
                j-=1

def generate():
    rcut = 12.0
    linecount = 0
    nspec=0
    with open('param.in', 'r') as inp:
        for line in inp:
            if linecount == 0:
                nspec = int(line.split()[0])
            linecount += 1
    param = param_init(nspec)
    linecount = 0
    with open('param.in', 'r') as inp:
        for line in inp:
            if linecount != 0:
               param.add_param(line.split()[0], line.split()[1], line.split()[2], line.split()[3],linecount-1)
            linecount += 1
    param.mix_param(nspec)
    param.print_genpot(nspec,rcut)
    
    
generate()

---------------------------------
---------------------------------
Generation of PSF Files:
---------------------------------
---------------------------------

CP2K needs a PSF topology file to generate the connectivity. These files are spacing sensitive and must be generated carefully.

These can be generated with VMD.

---------------------------------
---------------------------------
Constrained Dynamics: 
---------------------------------
---------------------------------


If you want to do a constraint, there are a couple of options: 

If you have something with three constraints - i.e. water:

g3x3 is the way to go - it works with three constraints (HH, OH, OH vectors).

The code block to add to the MOTION section is as follows:


 &CONSTRAINT
    CONSTRAINT_INIT T
    SHAKE_TOLERANCE 1E-10
    &G3X3
      DISTANCES 1.8897268 1.8897268 3.0859239
      MOLECULE 1
      ATOMS 1 2 3
    &END G3X3
  &END CONSTRAINT
  

  
---------------------------------
---------------------------------

Now if you want to do something else, like add a constraint between only two atoms you can use a collective variable. This requires adding a couple of other pieces of information.

First you need to add a COLVAR section to the SUBSYS section:


	&COLVAR
         &DISTANCE
             ATOMS 1 2
             &POINT
                 ATOMS 289
             &END POINT
             &POINT
                 ATOMS 290
             &END POINT
         &END DISTANCE
     &END COLVAR


Then you need to add a COLLECTIVE section in the CONSTRAINT section

    &COLLECTIVE
        COLVAR 1
        INTERMOLECULAR TRUE
        TARGET [angstrom] 2.6
    &END COLLECTIVE


And there you have it - a constraint between two separate atoms of separate molecules!

If you want to print out the value of the colvar - you need to do it via a sort of work around. You need to print it out through the metadynamics section of the code - but turning off the actually metadynamics i.e. DO_HILLS FALSE. The only thing that actually matter in this section are the METAVAR and the PRINT sections, but all keywords are needed.

  &FREE_ENERGY
    &METADYN
      DO_HILLS FALSE			! IMPORTANT
      DELTA_T 300.0
      NT_HILLS 400
      WELL_TEMPERED TRUE
      WW 0.001594
      &METAVAR					! Chooses which colvar to print
        COLVAR 1
        SCALE 1.0
      &END METAVAR
      &PRINT                   	! Prints the COLVAR out.
        &COLVAR
          &EACH
              MD 1
          &END EACH
          FILENAME =./out.colvar
        &END COLVAR
      &END PRINT
    &END METADYN
  &END FREE_ENERGY



---------------------------------
---------------------------------
The WHAM Method
---------------------------------
---------------------------------

I have also now developed a code that can set up a wham calculation. This is available in my code collection at https://github.com/piskuliche/My-Code-Collection.git

The trick is to use a harmonic restraint at different distances of about 0.01 hartree (but 0.005 in the code since it doesn't include the factor of 2 automatically.)

---------------------------------
---------------------------------
Ab Initio Molecular Dynamics:
---------------------------------
---------------------------------

The input file for AIMD differs mainly in the FORCE_EVAL section. I have included an example input file below. As you can see, the MOTION section is basically identical.

@SET data_path /people/cjmundy/SOURCE_FORGE/cp2k-6.1/cp2k/data
@SET BASE_NAME md
@SET ID 01
&GLOBAL
  PROJECT ${BASE_NAME}-${ID}
  RUN_TYPE MD
  PRINT_LEVEL LOW
  WALLTIME 23:50:00
&END GLOBAL
&MOTION
  &MD
    ENSEMBLE NVT
    STEPS 500000
    TIMESTEP 0.5
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE NOSE
      REGION MASSIVE
      &NOSE
        LENGTH              3
        YOSHIDA             3
        TIMECON  [wavenumber_t] 1000
        MTS                 2
      &END NOSE
    &END THERMOSTAT
    &PRINT
      &ENERGY
        FILENAME =${BASE_NAME}-${ID}.ener
        &EACH 
          MD 1
        &END EACH
      &END ENERGY   
    &END PRINT
  &END MD 
  &PRINT  
    &TRAJECTORY  SILENT
      &EACH 
        MD 1
      &END EACH
      FILENAME =${BASE_NAME}-${ID}.xyz
    &END TRAJECTORY  
    &VELOCITIES  SILENT
      &EACH 
        MD 1
      &END EACH
      FILENAME =${BASE_NAME}-${ID}.vel
    &END VELOCITIES 
    &RESTART
      FILENAME =restart
      &EACH 
        MD 1
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION

&FORCE_EVAL
  METHOD Quickstep                         				! GPW approach to AIMD
  &DFT
    WFN_RESTART_FILE_NAME ./RESTART						! Restart file to hold wfn ~1gb
    BASIS_SET_FILE_NAME ${data_path}/BASIS_MOLOPT		! Loc of Basis
    POTENTIAL_FILE_NAME ./POTENTIAL_revPBE				! Loc of Pot File
    &MGRID	
      CUTOFF 400										! Number of grid points - more = slower, but high accuracy
    &END MGRID
    &QS													! Quickstep parameters
      EPS_DEFAULT 1.0E-14								! basically a tolerance
      MAP_CONSISTENT									! Compute the exact derivative (Hks) of the energy with respect to the density matrix.
      EXTRAPOLATION ASPC								! Always Stable Predictor Corrector - MD Stability
      EXTRAPOLATION_ORDER 3								! Order of the extrapolation - higher = better, but more expensive and can affect stability
    &END QS
   &SCF   
     MAX_SCF            20								! Number of SCF iterations per step
     MAX_DIIS           7								! Number of vectors for extrapolation LC of errors
     EPS_SCF            1.0E-6							! Target accuracy.
     SCF_GUESS          RESTART							! Initial guess for w.f.
     &OUTER_SCF T    									! Predictor loop (with similar subparam meanings)
        EPS_SCF 1.0E-6			
        MAX_SCF 10  
     &END OUTER_SCF 
     &OT T												! Orbital Transformation section. 
       PRECONDITIONER   FULL_KINETIC					! Cholesky inversion of S and T, fast use for very large systems. 
       MINIMIZER        DIIS 							! Direct inversion in the iterative subspace: less reliable than CG, but sometimes about 50% faster
       N_DIIS           7								! Num of DIIS
     &END OT
     &PRINT SILENT										! turns off print
       &PROGRAM_RUN_INFO SILENT
       &END PROGRAM_RUN_INFO
     &END PRINT
   &END SCF												! Ends SCF
    &XC													! XC section
      &XC_GRID
        XC_SMOOTH_RHO  NN10                             ! How smoothing goes - method stuff (since grid - want continuity)
        XC_DERIV  SPLINE2_SMOOTH
      &END XC_GRID
      &XC_FUNCTIONAL 									! XC Functional
        &PBE
          PARAMETRIZATION REVPBE
        &END PBE
      &END XC_FUNCTIONAL
     &vdW_POTENTIAL
       DISPERSION_FUNCTIONAL PAIR_POTENTIAL
       &PAIR_POTENTIAL                                  ! Where the potential file is.
          TYPE DFTD3
          PARAMETER_FILE_NAME ${data_path}/dftd3.dat
          REFERENCE_FUNCTIONAL revPBE
       &END
     &END vdW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS												! This is the same.
    &CELL
      ABC 14.219066 14.219066 14.219066
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ./ion_pair_equil.xyz
    &END TOPOLOGY
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-revPBE-q6
    &END KIND
    &KIND Li
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-revPBE-q3
    &END KIND
    &KIND F
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-revPBE-q7
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-revPBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
!&EXT_RESTART
!  RESTART_FILE_NAME ./restart
!  RESTART_DEFAULT F
!  RESTART_POS T
!&END EXT_RESTART


---------------------------------
---------------------------------
Pseudopotential Optimization
---------------------------------
---------------------------------

To generate the necessary pseudopotentials, you need to run a different calculation. I have included an example input for Lithium, provided by chris.

&GLOBAL
  PROGRAM_NAME ATOM
&END GLOBAL
&ATOM
  ELEMENT Li										! The element in question

  RUN_TYPE PSEUDOPOTENTIAL_OPTIMIZATION				! Optimizes the pseudopotential

  ELECTRON_CONFIGURATION  [He] 2s1					
  CORE none											! Size of core electrons
  MAX_ANGULAR_MOMENTUM 2							! Largest ang mom calculated
  &METHOD		
     METHOD_TYPE  KOHN-SHAM							! DFT Method KS
     RELATIVISTIC DKH(2)							! Use Douglas-Kroll-Hess Hamiltonian of order 2 
     &XC
       &XC_FUNCTIONAL
          &PBE
             PARAMETRIZATION REVPBE					! Functional Name
          &END
       &END XC_FUNCTIONAL
     &END XC
  &END METHOD
  &OPTIMIZATION
    EPS_SCF 1.e-10									! Accuracy level
  &END

  &AE_BASIS
     BASIS_TYPE GEOMETRICAL_GTO						! type of basis all electron 
  &END AE_BASIS
  &PP_BASIS
     BASIS_TYPE GEOMETRICAL_GTO						! type of basis pseudopotential
  &END PP_BASIS
  &POTENTIAL
    PSEUDO_TYPE GTH									! pseudo type
    &GTH_POTENTIAL                                  ! Specifies the potential 
#Li GTH-PBE-q3 GTH-PBE
    3
     0.40000000    4   -14.08115455     9.62621962    -1.78361605     0.08515207
    0
    &END
    CONFINEMENT_TYPE  BARRIER						! Form of barrier
    CONFINEMENT 200. 4.0 12.0						! parameters.
  &END POTENTIAL

  &POWELL											! Powell optimization
     ACCURACY   1.e-14
     STEP_SIZE  0.08
     MAX_INIT   5
     MAX_FUN    1000
     STEP_SIZE_SCALING 0.80
     WEIGHT_PSIR0 0.0
     TARGET_POT_SEMICORE      [eV]      0.003000
     TARGET_POT_VALENCE       [eV]      0.000100
     TARGET_POT_VIRTUAL       [eV]      0.00100
     WEIGHT_POT_NODE                    10.0
     WEIGHT_POT_SEMICORE                2.0
     WEIGHT_POT_VALENCE                 5.0
     WEIGHT_POT_VIRTUAL                 1.0
     SEMICORE_LEVEL           [eV]      20.0
  &END

&END ATOM


