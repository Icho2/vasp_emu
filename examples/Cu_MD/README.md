# Example: Molecular Dynamics on a Copper System

This example focuses on relaxing the bond length of an oxygen dimer through geometry optimization with the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm using a potential from Open Catalyst Project (OCP).

## Input

Only two files are needed to run this example:

### POSCAR

The POSCAR file contains out system of interest, crystal composed of 108 Copper atoms. It follows the same formatting as you would see in VASP.

### INCAR

```Text
EDIFFG = -0.01
ISIF = 1
MD_ALGO = 0
NSW = 1000
IBRION=0
IOPT=0
MAXMOVE=0.1
POTENTIAL=EMT
POTIM=5.0
```

The INCAR file contains settings for our emulator. In this example, `IOPT` and `POTENTIAL` are the only tags that are not part of VASP. `IOPT=0` tells the emulator that `IBRION` specifies the optimizer, like in VASP. In this example, `IOPT=0` and `IBRION=0`. This means that this example will perform Molecular Dynamics, specifically using the Velocity-Verlet algorithm due to `MD_ALGO=0`. `POTENTIAL` tells the emulator which potential to use. In this example, `POTENTIAL=EMT`for simplicity. 

## Calculation and Expected Output


To run this example, you only need to run `vasp_emu`. Several files will be generated:

- CONTCAR: the final structure
- OSZICAR: a copy of all the information printed to standard out
- OUTCAR: contains job-specific information
- XDATCAR: contains all the structures during the optimization
- dynamics.log: produced by the optimizer
- md.traj: Trajectory file that shows the molecular simulation

