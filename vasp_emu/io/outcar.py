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

    
    def write_header(self, atoms: Atoms, neb_params: dict = None) -> None:
        """
        Write the header for the OUTCAR file.
        Conditionally writes NEB header info if provided.
        """
        self.info("                                       "
                  f"number of ions     NIONS ={len(atoms):7d}\n")
        
        # If this is an NEB job, the NEBJob class will pass neb_params
        if neb_params:
            self.write_neb_header(**neb_params)


    def write_neb_header(self, spring: float, lclimb: bool, efirst: float, elast: float, **kwargs) -> None:
        """
        Writes the NEB-specific header block, formatted to match VASP's output.
        """
        # These values are not in config.yaml, so we use their VASP defaults from neb.F
        ltangentold = kwargs.get('ltangentold', False)
        ldneb = kwargs.get('ldneb', False)
        ldneborg = kwargs.get('ldneborg', False)
        lnebcell = kwargs.get('lnebcell', False)

        # Helper for Fortran-style boolean formatting ('T' or 'F')
        def format_bool(b):
            return 'T' if b else 'F'

        self.info(f" CHAIN: Read ICHAIN         0") # As in VASP
        self.info(f" CHAIN: Running the NEB")
        self.info(f" NEB:        SPRING    {spring:14.6f}")
        self.info(f" NEB:        LCLIMB       {format_bool(lclimb):>7}")
        self.info(f" NEB:   LTANGENTOLD       {format_bool(ltangentold):>7}")
        self.info(f" NEB:         LDNEB       {format_bool(ldneb):>7}")
        self.info(f" NEB:      LDNEBORG       {format_bool(ldneborg):>7}")
        self.info(f" NEB:      LNEBCELL       {format_bool(lnebcell):>7}")
        self.info(f" NEB:        EFIRST    {efirst:14.6f}")
        self.info(f" NEB:         ELAST    {elast:14.6f}")


    def write_step(self, step: int, step_data: dict, neb_stats: dict = None) -> None:
        """
        Writes a full ionic step to the OUTCAR using only pre-calculated data.
        
        Args:
            step (int): The current ionic step number.
            step_data (dict): A dictionary from Job.get_step_data().
            neb_stats (dict, optional): Dict with NEB-specific stats.
        """
        self.write_ionic_header(step)

        if neb_stats:
            self.write_neb_step_stats(**neb_stats)

        self.write_common_stats(step_data) # Pass the whole dict for convenience


    def write_ionic_header(self, step: int) -> None:
        """
        Write the header for an ionic step in the OUTCAR file.

        Args:
            step (int): The current ionic step number.
        """
        self.info(f'\n--------------------------------------- Ionic step {step:>8}  -------------------------------------------\n\n\n\n')


    def write_neb_step_stats(self, dist_prev: float, dist_next: float, angle: float, 
                             proj_spring: float, proj_real: float) -> None:
        """
        Writes the NEB-specific information for an ionic step.
        """
        line1 = (f"  NEB: distance to prev, next image, angle between "
                 f" {dist_prev:11.6f} {dist_next:11.6f} {angle:11.6f}")  # is it 11.6f or 12.6f?
        line2 = (f"  NEB: projections on to tangent (spring, REAL) "
                 f" {proj_spring:11.6f} {proj_real:11.6f}")  # is it 11.6f or 12.6f?
        self.info(line1)
        self.info(line2)


    def write_common_stats(self, data: dict) -> None:
        """
        Writes the common statistics block from a data dictionary.
        """
        self.write_pos_forces(data['positions'], data['forces'])
        self.write_energy(data['energy'])
        self.info(f"  FORCES: max atom, RMS {data['fmax_atom']:11.6f} {data['f_rms']:11.6f}")
        self.info(f"  Stress total and by dimension {data['stress_total']:11.6f} {data['stress_dim']:11.6f}")
        self.info(f"  volume of cell : {data['volume']:11.2f}")
        self.info(f" number of electron {data['n_electrons']:15.7f} magnetization {data['magnetization']:15.7f}")


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

        # NOTE: yes, there are 2 spaces between 'energy' and 'without', and none between 'entropy' and '='
        line = f"  energy  without entropy= {energy:17.8f}  energy(sigma->0) = {energy:17.8f}\n"
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