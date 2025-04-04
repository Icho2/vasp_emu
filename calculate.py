from vasp import read_incar
import ase
import numpy as np
import os
from ase.calculators.emt import EMT
import sys
from tsase.optimize import SDLBFGS

#Read incar
incar = read_incar()
iopt = incar['iopt']
ibrion = incar['ibrion']

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
poscar.calc = potential

#set up optimizer 
calc = None

if iopt == 1 and ibrion != 1:
    calc = SDLBFGS(poscar, trajectory='path.traj', logfile='optimization.log', maxstep=incar['maxstep'], damping=incar['damping'])

elif iopt == 7:
    from ase.optimize import FIRE
    calc = FIRE(poscar, trajectory='path.traj', logfile='optimization.log', maxstep=incar['maxstep'], dtmax=incar['dtmax'], Nmin=incar['nmin'], finc=incar['finc'], fdec=incar['fdec'], astart=incar['astart'], fa=incar['fa'], a=incar['a'])

elif ibrion == 1 and iopt != 1:
    from ase.optimize import BFGS 
    calc = BFGS(poscar, trajectory='path.traj', logfile='optimization.log', maxstep=incar['maxstep'], alpha=incar['alpha'])

elif ibrion == 2:
    sys.exit("dont waste your time here")

elif iopt == 3 and ibrion == 3:
    from ase.optimize import MDMin
    calc = MDMin(poscar, trajectory='path.traj', logfile='optimization.log', dt=incar['dt'], maxstep=incar['maxstep'])

elif ibrion == 1 and iopt == 1:
    sys.exit('please change your iopt or ibrion')
else:
    sys.exit('please set a calculation')

max_steps=incar['nsw']
steps=0
converged=False
outcar = open("OUTCAR", 'w')
while not converged:
        calc.run(fmax=-1*incar['ediffg'],steps=1)
        outcar.write(f'U: {poscar.get_potential_energy()}\n')
        steps+=1
        if steps==max_steps:
           sys.exit('we reached NSW')
        if get_fmax(poscar)<-1*incar['ediffg']:
           converged = True        	

outcar.close()
# Making XDATCAR and CONTCAR files
traj = ase.io.trajectory.Trajectory('path.traj');ase.io.vasp.write_vasp_xdatcar("XDATCAR", traj);os.remove('path.traj')
ase.io.write('CONTCAR',poscar,format='vasp')
