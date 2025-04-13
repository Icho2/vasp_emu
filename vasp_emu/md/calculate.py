from ase.md import MDLogger
from ase.io import read
from ase import units

def calculate(incar):
    from ase.md.verlet import VelocityVerlet
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    atoms = read('POSCAR')
    
    if incar['potential'] == 'OCP':
        from fairchem.core import OCPCalculator
        potential = OCPCalculator(
            model_name='EquiformerV2-31M-S2EF-OC20-All+MD',
            local_cache="pretrained_models",
            seed=1234,
            cpu=True,
        )
    if incar['potential'] == 'LJ':
        from ase.calculators.lj import LennardJones
        potential = LennardJones()

    if incar['potential'] == 'EMT':
        from ase.calculators.emt import EMT
        potential = EMT()

    atoms.calc = potential

    MaxwellBoltzmannDistribution(atoms, temperature_K=300)
    def printenergy(a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

    from ase.md.verlet import VelocityVerlet
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

    dyn = VelocityVerlet(atoms=atoms, timestep=incar['potim'], trajectory='MD.traj', logfile='optimization.log')
    dyn.attach(printenergy, interval=10)
    dyn.run(incar['nsw'])
