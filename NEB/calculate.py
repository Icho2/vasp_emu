from ase.io import read, write
from ase.mep import NEB
from ase.optimize import MDMin
from ase import Atoms
from fairchem.core import OCPCalculator
from vasp import read_incar

def calculate(incar):
    # Read intial and final states:
    initial = read('initial.traj')
    final = read('final.traj')

    # Make a band consiting of 5 images:
    images = [initial]
    images += [initial.copy() for i in range(3)]
    images += [final]

    # Sort atoms by symbols and positions
    for image in images:
        image.arrays['positions'] = image.arrays['positions'][image.arrays['numbers'].argsort()]
        image.arrays['numbers'] = image.arrays['numbers'][image.arrays['numbers'].argsort()]

    for i, atoms in enumerate(images):
        print(f"Image {i}:")
        for j, atom in enumerate(atoms):
            print(f"Atom {j}: {atom.symbol}, Position: {atom.position}")

    neb = NEB(images, allow_shared_calculator=True)

    # Interpolate linearly the positions of the three middle images:
    neb.interpolate()

    if incar['potential'] == 'OCP':
        potential = OCPCalculator(
            model_name="EquiformerV2-31M-S2EF-OC20-All+MD",
            local_cache="pretrained_models",
            cpu=True,
        )

    elif incar['POTENTIAL'] == 'Vasp':
        potential = Vasp()

    # Set calculators:
    for image in images[1:4]:
        image.calc = potential

    # Optimize
    optimizer = MDMin(neb, trajectory='A2B.traj')
    optimizer.run(fmax=0.04)
