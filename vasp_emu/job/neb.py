"""Modules"""
import sys
import time
from math import sqrt

from ase.mep import NEB
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

class NEBJob(Job):
    """ 
    An instance of the Job class used to perform nudged elastic band method
    
    Attributes:
        job_name (str): name of the job
        images (list of Atoms) : images that make up the band
        neb (NEB) : NEB object to drive the dynamics
    """
    def __init__(self,**kwargs):
        """
        Construct an NEBJob
        
        Arguments:
            Please refer to the parent Job Class
        """
        super().__init__(**kwargs)
        self.job_name = "nudged-elastic-band"
        self.images = [self.structures["initial"]]
        self.images.extend([self.structures["initial"].copy()]*self.job_params['num_img'])
        self.images.extend([self.structures["final"]])
        # Sort atoms by symbols and positions
        for image in self.images:
            image.arrays['positions'] = image.arrays['positions'][image.arrays['numbers'].argsort()]
            image.arrays['numbers'] = image.arrays['numbers'][image.arrays['numbers'].argsort()]
        if self.logger is not None:
            for i, atoms in enumerate(self.images):
                self.logger.info(f"Image {i}:")
                for j, atom in enumerate(atoms):
                    self.logger.info(f"Atom {j}: {atom.symbol}, Position: {atom.position}")
        self.neb = NEB(self.images, allow_shared_calculator=True)
        self.neb.interpolate()


    def set_dynamics(self, name) -> None:
        """
        Extend the set_dynamics function in the parent Job class by modifying the logger
        Redirect output to both stdout and a file
        
        Arguments:
            name (str): the input argument for set_dynamics in the parent Job class
        """
        super().set_dynamics(name)
        self.dynamics.log = opt_log.__get__(self.dynamics,Optimizer)
        self.dynamics.attach(lambda : self.dyn_logger.info(self.dynamics.log()),interval=1)

    def calculate(self) -> None:
        """
        Perform the NEB Calculation
        """
        for image in self.images:
            image.calc = self.potential
        max_force = self.job_params["fmax"]
        max_steps = self.job_params["max_steps"]
        
        self.dynamics.atoms = self.neb
        steps = 0
        finished = False
        while not finished:
            self.dynamics.run(fmax=max_force,steps=1)
            steps += 1
            if steps == max_steps:
                if self.logger is not None:
                    self.logger.info('Reached NSW')
                finished = True
            # if self.get_fmax(self.curr_structure)<fmax:
                # finished = True
        # self.create_xdatcar()
