"""Modules"""
import time
from math import sqrt

import ase
from ase.constraints import FixAtoms
from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate

from vasp_emu.job.job import Job


def opt_log(self, forces=None):
    """
        Redefine the behavior of the log function for Optimizers. 
        
        Arguments:
            forces (list of floats) : The forces that correspond to the optimizable object
        Returns:
            Message to be sent to the logger
    """
    if forces is None:
        forces = self.optimizable.get_forces()
    fmax = sqrt((forces ** 2).sum(axis=1).max())
    e = self.optimizable.get_potential_energy()
    T = time.localtime()
    name = self.__class__.__name__
    # everything above this line exactly matches the Optimizer.log()
    msg = ""
    if self.nsteps == 0:
        msg += "=======================================================\n"
        msg += f"{" " * len(name)}  {"Step":4s} {"Time":>9s} {"Energy":>13s}  {"fmax":>10s}\n"

    msg += f"{name}:  {self.nsteps:3d}    {t[3]:02d}:{t[4]:02d}:{t[5]:02d} {e:12.6f} {fmax:12.6f}"
    return msg

class DimerJob(Job):
    """ 
        An instance of the Job class used to run dimer
        
        Attributes:
            job_name (str): name of the job
    """
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.job_name = "dimer"

    def calculate(self) -> None:
        """
        Perform the dimer calculation
        """
        curr_structure = self.structures["initial"]
        mask = [atom.tag > 0 for atom in curr_structure]
        curr_structure.set_constraint(FixAtoms(mask=mask))
        curr_structure.calc = self.potential
        curr_structure.get_potential_energy() #为什么?
        max_force = self.job_params["fmax"]
        max_steps = self.job_params["max_steps"]

        steps = 0
        finished = False

        with DimerControl(initial_eigenmode_method='displacement', displacement_method='vector', logfile=None, mask=[0,0,0,0,1]) as d_control:
            d_atoms = MinModeAtoms(curr_structure,d_control)

            # Displace the atoms
            displacement_vector = [[0.0] * 3] * 5
            displacement_vector[-1][1] = -0.1
            d_atoms.displace(displacement_vector=displacement_vector)

            with MinModeTranslate(d_atoms, trajectory='dynamics.traj') as dim_rlx:
                while not finished:
                    dim_rlx.run(fmax=max_force,steps=1)
                    if self.logger is not None:
                        self.logger.debug(f'U: {curr_structure.get_potential_energy()}' + \
                                            f'   fmax: {self.get_fmax(curr_structure)}')
                    # CONTCAR should be written after each step, used to restart jobs
                    ase.io.write('CONTCAR',curr_structure,format='vasp')
                    steps+=1
                    if self.get_fmax(curr_structure) < max_force:
                        finished = True
                    elif steps == max_steps:
                        if self.logger is not None:
                            self.logger.info('Reached NSW')
                        finished  = True
        # self.create_xdatcar()
