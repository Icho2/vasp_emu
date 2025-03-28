import numpy
import ase
from ase.build import molecule
import tsase

#a = ase.io.read("al.traj")
a = molecule('H2')
al = tsase.calculators.lepsho()
a.set_calculator(lepsho)

u = a.get_potential_energy()
print(u, "eV")

f = a.get_forces()
print(numpy.linalg.norm(f), "eV/angstrom")
