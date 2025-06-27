"""Modules"""
import time
from math import sqrt

import ase
from ase.optimize.optimize import Optimizer

from vasp_emu.job.job import Job

def opt_log(self, forces=None) -> str:
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
    t = time.localtime()
    name = self.__class__.__name__
    # everything above this line exactly matches the Optimizer.log()
    msg = ""
    if self.nsteps == 0:
        msg += "=======================================================\n"
        msg += f'{" " * len(name)}  {"Step":4s} {"Time":>9s} {"Energy":>13s}  {"fmax":>10s}\n'

    msg += f"{name}:  {self.nsteps:3d}    {t[3]:02d}:{t[4]:02d}:{t[5]:02d} {e:12.6f} {fmax:12.6f}"
    return msg

class OptJob(Job):
    """ 
    An instance of the Job class used to perform geometry optimization
    
    Attributes:
        job_name (str): name of the job
    """
    def __init__(self,**kwargs):
        """
        Construct an OptJob
        
        Arguments:
            Please refer to the parent Job Class
        """
        super().__init__(**kwargs) # always first
        self.job_name = "optimization"
        self.set_dynamics() # always last

    def set_dynamics(self) -> None:
        """
        Extend the set_dynamics function in the parent Job class by modifying the logger
        Redirect output to both stdout and a file
        """
        super().set_dynamics()
        self.dynamics.log = opt_log.__get__(self.dynamics,Optimizer)
        self.dynamics.attach(lambda : self.dyn_logger.info(self.dynamics.log()),interval=1)

    def calculate(self) -> None:
        """
        Perform the geometry optimization
        """
        curr_structure = self.poscar
        max_force = self.job_params["fmax"]
        max_steps = self.job_params["max_steps"]
        curr_structure.calc = self.potential
        steps = 0
        finished = False
        while not finished:
            self.dynamics.run(fmax=max_force,steps=1)
            self.logger.info(f'U: {curr_structure.get_potential_energy()}   ' + \
                                f'fmax: {self.get_fmax(curr_structure)}')
            # CONTCAR should be written after each step, used to restart jobs
            ase.io.write('CONTCAR',curr_structure,format='vasp')
            steps += 1
            if self.get_fmax(curr_structure) < max_force:
                finished = True
            elif steps == max_steps:
                self.logger.info('Reached NSW')
                finished = True
        self.create_xdatcar()
