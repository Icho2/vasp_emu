import numpy as np

class double_well():
    def __init__(self, a=4.5, b=3.5):
        self.a = a
        self.b = b
        self.atoms = None
        self.u = None
        self.f = None

    def calculate(self):
        pos = self.atoms.positions
        forces = pos * 0.0
        forces[0] = (self.a*4*(pos[0]**4))- (3*(pos[0]**2)) - (2*self.b*pos[0])
        forces[1] = 0.0
        forces[2] = 0.0
        self.f = forces
        self.u = (self.a*(pos[0]**4)) - (pos[0]**3) - (self.b*(pos[0]**2))

    def get_potential_energy(self, atoms, force_consistent=False):
        if self.calculation_required(atom, "energy"):
            self.atoms = atoms.copy()
            self.calculate()
        return self.u
    
    def get_forces(self, atoms):
        if self.calculation_required(atoms, "forces"):
            self.atoms = atoms.copy()
            self.calculate()
        return self.f.copy()

    def get_stress(self, atoms):
        return np.zeros((1))
    
    def calculation_required(self, atoms, quantities):
        if atoms != self.atoms or self.atoms == None:
            return True
        if np.all(self.f) == None or np.all(self.u) == None or np.all(atoms) == None:
            return True
        return False

    def set_atoms(self, atoms):
        pass
