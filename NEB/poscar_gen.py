from ase.io import write, read
from ase.build import molecule
from ase import Atoms

# Make the atom objects
CH3OH = molecule('CH3OH')
CH3Cl = molecule('CH3Cl')
NaOH = Atoms('NaOH', positions=[[1.0,1.0,1.0],[1.0,2.0,1.0],[1.0,2.0,2.0]])
NaCl = molecule('NaCl')
reactants = CH3Cl + NaOH
products = CH3OH + NaCl

write('initial.traj', reactants)
write('final.traj', products)

initial = read('initial.traj', index=':')
final = read('final.traj', index=':')
initial.sort()
final.sort()

for i, atoms in enumerate(initial):
    print(f"Image {i}:")
    for j, atom in enumerate(atoms):
        print(f"Atom {j}: {atom.symbol}, Position: {atom.position}")

for i, atoms in enumerate(final):
    print(f"Image {i}:")
    for j, atom in enumerate(atoms):
        print(f"Atom {j}: {atom.symbol}, Position: {atom.position}")

write('initial.traj', initial)
write('final.traj', final)
