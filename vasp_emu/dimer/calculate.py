import pytest

from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import read
from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate

def test_dimer_method(testdir):
    atoms = read('POSCAR')
    
    # Freeze System
    mask = [atom.tag > 0 for atom in atoms]
    atoms.set_constraint(FixAtoms(mask=mask))

    # Use calculator
    atoms.calc = EMT()
    atoms.get_potential_energy()
    
    # Setting up Dimer
    with DimerControl(initial_eigenmode_method='displacement', displacement_method='vector', logfile=None, mask=[0,0,0,0,1]) as d_control:
        d_atoms = MinModeAtoms(atoms, d_control)
        
        # Displace the atoms
        displacement_vector = [[0.0] * 3] * 5
        displacement_vector[-1][1] = -0.1
        d_atoms.displace(displacement_vector=displacement_vector)
        
        # Converge to a saddle point 
        with MinModeTranslate(d_atoms, trajectory='dimer_method.traj', logfile=None) as dim_rlx:
            dim_rlx.run(fmax=0.001)
