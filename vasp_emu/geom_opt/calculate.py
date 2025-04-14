from vasp_emu.vasp import read_incar
import ase
import numpy as np
import os
import sys
from tsase.optimize import SDLBFGS

# Vasp
from ase.calculators.vasp import Vasp

def calculate(incar):
    # Reading INCAR tags
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

    if incar['potential'] == 'OCP':
        from fairchem.core import OCPCalculator
        potential = OCPCalculator(
            model_name='EquiformerV2-31M-S2EF-OC20-All+MD', 
            local_cache="pretrained_models",
            seed=1234,
            cpu=True,
        )
    
    elif incar['potential'] == 'VASP':
        from ase.calculators.vasp import Vasp
        potential = Vasp(atoms=poscar,
            command='mpirun vasp_std',
        )
        
    elif incar['potential'] == 'PYAMFF':
        from pyamff.ase_calc_fortran import aseCalcF
        try:
            f = open('mlff.pyamff')
            potential = aseCalcF()
        except:
            from pyamff.ase_calc import aseCalc
            try:
                f = open('pyamff.pt')
                potential = aseCalc('pyamff.pt')
            except:
                sys.exit('Please add a model.')
        
        f.close()
        
    poscar.calc = potential

    # Set up optimizer 
    dyn = None

    if iopt == 1 and ibrion != 1:
        dyn = SDLBFGS(poscar, trajectory='path.traj', logfile='optimization.log', maxstep=incar['maxstep'], damping=incar['damping'])

    elif iopt == 7:
        from ase.optimize import FIRE
        dyn = FIRE(poscar, trajectory='path.traj', logfile='optimization.log', maxstep=incar['maxstep'], dtmax=incar['dtmax'], Nmin=incar['nmin'], finc=incar['finc'], fdec=incar['fdec'], astart=incar['astart'], fa=incar['fa'], a=incar['a'])

    elif ibrion == 1 and iopt != 1:
        from ase.optimize import BFGS 
        dyn = BFGS(poscar, trajectory='path.traj', logfile='optimization.log', maxstep=incar['maxstep'], alpha=incar['alpha'])

    elif ibrion == 2:
        sys.exit("dont waste your time here")

    elif iopt == 3 and ibrion == 3:
        from ase.optimize import MDMin
        dyn = MDMin(poscar, trajectory='path.traj', logfile='optimization.log', dt=incar['dt'], maxstep=incar['maxstep'])

    elif ibrion == 1 and iopt == 1:
        sys.exit('please change your iopt or ibrion')
    else:
        sys.exit('please set a calculation')

    max_steps=incar['nsw']
    steps=0
    converged=False
    outcar = open("OUTCAR", 'w')
    while not converged:
            dyn.run(fmax=-1*incar['ediffg'],steps=1)
            outcar.write(f'U: {poscar.get_potential_energy()}   fmax: {get_fmax(poscar)}\n')
            steps+=1
            if steps==max_steps:
                sys.exit('we reached NSW')
            if get_fmax(poscar)<-1*incar['ediffg']:
                converged = True        	

    outcar.close()
    
    # Making XDATCAR and CONTCAR files
    traj = ase.io.trajectory.Trajectory('path.traj');ase.io.vasp.write_vasp_xdatcar("XDATCAR", traj);os.remove('path.traj')
    ase.io.write('CONTCAR',poscar,format='vasp')
    
    
    return

def get_fmax(poscar, **kwargs):
        """Returns fmax, as used by optimizers with NEB."""
        forces = poscar.get_forces()
        return np.sqrt((forces ** 2).sum(axis=1).max())

