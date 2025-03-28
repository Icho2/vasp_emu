from vasp import read_incar
import ase
import numpy as np
from ase.calculators.emt import EMT
import sys
#from tsase.tsase.optimize import SDLBFGS
from tsase.optimize import SDLBFGS

#Read incar
incar = read_incar()
iopt = incar['iopt']
ibrion=incar['ibrion']
print(incar)

#Read poscar
poscar = ase.io.read('POSCAR')

#calculate with ase calculator
'''
1. Set optimizer
2. Set potential (calculator)
3. Convergence criteria
4. Max steps
'''

def get_fmax(poscar, **kwargs):
        """Returns fmax, as used by optimizers with NEB."""
        forces = poscar.get_forces()
        return np.sqrt((forces ** 2).sum(axis=1).max())

potential = EMT()
#print(incar['nsw'])
poscar.calc = potential

#set up optimizer 
calc = None

if iopt == 1:
    calc = SDLBFGS(poscar)

elif iopt == 7:
    from ase.optimize import FIRE
    calc = FIRE(poscar)

elif ibrion == 1:
    from ase.optimize import BFGS 
    calc = BFGS(poscar)

elif ibrion == 2:
    sys.exit("dont waste your time here")
    #from tsase import CG 
    #opt = cg

elif ibrion == 3:
    from ase.optimize import MDMin
    calc = MDMin(poscar)

else:
    sys.exit('please set a calculation')

#calc = BFGS(poscar)
max_steps=incar['nsw']
steps=0
converged=False
outcar = open("OUTCAR", 'w')
while not converged:
        calc.run(fmax=-1*incar['ediffg'],steps=1)
        #if fmax<-1*incar['ediffg']:
        #   converged=True
        outcar.write(f'U: {poscar.get_potential_energy()}\n')
        steps+=1
        if steps==max_steps:
           sys.exit('we reached NSW')
        if get_fmax(poscar)<-1*incar['ediffg']:
           converged = True        	

outcar.close()
ase.io.write('CONTCAR',poscar,format='vasp')

