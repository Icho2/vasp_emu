# Example: O2 Relaxation

This example focuses on relaxing the bond length of an oxygen dimer through geometry optimization with the Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm using a PyAMFF. In this folder, you will notice a directory called `from_pyamff`. This folder contains the config.ini and trajectory file (train.traj) used to train the PyAMFF potential. It also contains the final potential files produced by PyAMFF: mlff.pyamff and pyamff.pt. This folder is not required to run the emulator. It is merely listed here as a reference for any user to gain information about the training parameters that were used.

## Input

Only three files are needed to run this example:

 - ### POSCAR

```Text
H  Pd H 
 1.0000000000000000
    40.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   40.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   40.0000000000000000
 H   Pd  H  
   1   3   1
Cartesian
 20.9528395179931124 20.1758617195403609 18.5939685206410665
 21.0363745179931101 19.4146347195403628 20.0546385206410669
 19.9302309830626143 20.5805072330725984 21.4060314793589335
 18.9636254820068899 18.9414275902822737 20.0011362814718474
 19.9335445748133750 21.0585724097177263 20.0493856662234933
```

The POSCAR file contains out system of interest. It follows the same formatting as you would see in VASP.

### INCAR

```Text
EDIFFG = -0.01
POTIM = 1
ISIF = 1
MD_ALGO = 0
NSW = 100
IBRION=3
IOPT=1
MAXMOVE=0.1
POTENTIAL=PYAMFF
```

The INCAR file contains settings for our emulator. In this example, `IOPT` and `POTENTIAL` are the only tags that are not part of VASP. `IOPT=1` tells the emulator to use the SDLBFGS optimizer. NOTE that in VASP, the "RMM-DIIS" optimizer would be used for the same settings. `POTENTIAL` tells the emulator which machine-learned potential to use. In this example, `POTENTIAL=PYAMFF`, so a potential from PYAMFF is being used.

- ### mlff.pyamff or pyamff.pt

If you are familiar with PyAMFF, then you know that the PyAMFF potentials are written to a binary file: pyamff.pt as well as a text file: pyamff.pt. Either one of these files can be used to load the potential, and in fact, the emulator will use the first one that it finds.

## Calculation and Expected Output

During the calculation, a temporary file called "path.traj" will be generated. Therefore, if you already have a file with this name in your directory, it's recommended that you temporarily rename.

To run this example, you only need to run `vasp_emu`. Several files will be generated:

- CONTCAR: the final structure
- OSZICAR: a copy of all the information printed to standard out
- OUTCAR: contains job-specific information
- XDATCAR: contains all the structures during the optimization
- dynamics.log: produced by the optimizer

NOTE: Since we are using an OCP Potential, there will also be a folder generated called "pretrained_models". This folder holds the OCP Potential.
