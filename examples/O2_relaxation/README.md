# Example: O2 Relaxation

This example focuses on relaxing the bond length of an oxygen dimer through geometry optimization with the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.

## Input

Only two files are needed to run this example:

### POSCAR

```Text
O 
 1.0000000000000000
     8.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    8.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    8.0000000000000000
 O  
   2
Cartesian
  4.0000000000000000  4.0000000000000000  3.3827007447170030
  4.0000000000000000  4.0000000000000000  5.8172992552829976
```

The POSCAR file contains out system of interest. It follows the same formatting as you would see in VASP.

### INCAR

```Text
EDIFFG = -0.01
POTIM = 1
ISIF = 1
MD_ALGO = 0
NSW = 100
IBRION=1
IOPT=0
MAXMOVE=0.1
POTENTIAL=OCP
```

The INCAR file contains settings for our emulator. In this example, `IOPT` and `POTENTIAL` are the only tags that are not part of VASP. `IOPT` tells the emulator to `IBRION` specifies the optimizer, like in VASP. In this example, `IOPT=0` and `IBRION=1`. This means that this example will use the BFGS example. NOTE that in VASP, the "RMM-DIIS" optimizer would be used for the same settings. `POTENTIAL` tells the emulator which machine-learned potential to use. In this example, `POTENTIAL=OCP`, so a potential from the Open Catalyst Project (OCP) is being used. Since the INCAR does not also specify which potential, the default ("EquiformerV2-31M-S2EF-OC20-All+MD") is being used.

## Calculation and Expected Output

During the calculation, a temporary file called "path.traj" will be generated. Therefore, if you already have a file with this name in your directory, it's recommended that you temporarily rename.

To run this example, you only need to run `vasp_emu`. Several files will be generated:

- CONTCAR: the final structure
- OSZICAR: a copy of all the information printed to standard out
- OUTCAR: contains job-specific information
- XDATCAR: contains all the structures during the optimization
- optimization.log: produced by the optimizer

NOTE: Since we are using an OCP Potential, there will also be a folder generated called "pretrained_models". This folder holds the OCP Potential.