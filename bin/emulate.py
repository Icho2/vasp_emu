#!/usr/bin/env python3
# Here we will read the INCAR file and decide what the user wants to run. It will be the one and only executable all users run.
import subprocess

from vasp import read_incar

#read incar
incar = read_incar()


if incar['emulate'] == 0: # geometry_relaxation
    from geom_opt.calculate import calculate
    calculate(incar)
    
elif incar['emulate'] == 1: #NEB
    from NEB.calculate import calculate
    calculate(incar)

elif incar['emulate'] == 2: #Dimer
    from DIMER.calculate import test_dimer_method
    test_dimer_method(incar)
    
else:
    print('Fix emulate tag')
