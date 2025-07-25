"""Handles formatting and writing of OUTCAR"""
from ase.atoms import Atoms
from logging import Logger, FileHandler, Formatter
import numpy as np


class OutcarWriter(Logger):
    """
    A class to handle writing to the OUTCAR file.
    Inherits from Logger to use logging functionality.
    """
    def __init__(self, dir: str = '.') -> None:
        """
        Initialize the OutcarWriter at `dir/outcar`.

        Args:
            dir (str): Directory where the OUTCAR file will be written.
        """
        super().__init__(name=dir + '/outcar')
        # Prevent adding multiple handlers if the logger instance is somehow reused by name
        if not self.handlers:
            self.propagate = False
            self.setLevel('INFO')
            file_handler = FileHandler(filename=f"{dir}/OUTCAR", mode='w')
            file_handler.setFormatter(Formatter('%(message)s'))
            self.addHandler(file_handler)


    def write(self, atoms:Atoms, step: int) -> None:
        """
        A method that calls all other write methods below
        
        Args:
            atoms (Atoms): The Atoms object to write to OUTCAR.
            step (int): The current ionic step number.
        """
        self.write_ionic_step_header(step)
        self.write_pos_forces(atoms.get_positions(), atoms.get_forces())
        self.write_energy(atoms.get_potential_energy())


    def write_ionic_step_header(self, step: int) -> None:
        """
        Write the header for an ionic step in the OUTCAR file.

        Args:
            step (int): The current ionic step number.
        """
        self.info(f'\n--------------------------------------- Ionic step {step:>8}  -------------------------------------------\n\n\n\n')


    def write_energy(self, energy: float) -> None:
        """
        Write the energy to the OUTCAR file.
        
        Args:
            energy (float): The energy to write.
        """
        if energy <= -1e6 or energy >= 1e7:
            self.warning(f"Energy {energy} is out of expected range for vtstscripts\n"
                         "(expected between characters 67 ad 78)")
        if energy <= -1e9 or energy >= 1e10:
            self.warning(f"Energy {energy} is out of expected range for proper printing in OUTCAR")
        line = f"  energy without entropy ={energy:18.8f}  energy(sigma->0) ={energy:18.8f}\n"
        self.info(line)


    def write_pos_forces(self, pos: np.ndarray, forces: np.ndarray) -> None:
        """
        Write the position and forces block to the OUTCAR file.
        
        Args:
            pos (np.ndarray): The positions of atoms (N x 3).
            forces (np.ndarray): The forces on atoms (N x 3).
        """
        self.info(" POSITION                                       TOTAL-FORCE (eV/Angst)")
        self.info(" -----------------------------------------------------------------------------------")
        assert pos.shape == forces.shape, "Position and forces must have the same shape"
        for i, (p, f) in enumerate(zip(pos, forces)):
            line = f"{p[0]:13.5f} {p[1]:13.5f} {p[2]:13.5f}  {f[0]:14.6f} {f[1]:14.6f} {f[2]:14.6f}"  # XXX: do we need to subtract drift?
            self.info(line)
        self.info(" -----------------------------------------------------------------------------------\n")