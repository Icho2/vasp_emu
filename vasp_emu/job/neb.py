"""Modules"""
import sys, os
import time
from math import sqrt

import ase
from ase.mep import NEB
from ase.optimize.optimize import Optimizer, OptimizableAtoms
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

        # make band
        n_images = self.job_params['num_img']
        self.images = []
        interpolate = False if self.structures["final"] is None else True
        if interpolate:
            # final image was provided; construct intermediates now
            self.images = [self.structures["initial"]]
            for _ in range(n_images):
                self.images.append(self.structures["initial"].copy())
            self.images.append(self.structures["final"])
        else:
            # assume nebmake was already run; 00, 01, ... directories exist
            for i in range(n_images+2):
                i_dir = '0'+str(i) if i < 10 else str(i)
                poscar = os.path.join(i_dir, "POSCAR")
                curr_structure = ase.io.read(poscar)
                self.images.append(curr_structure)

        if self.logger is not None:
            for i, atoms in enumerate(self.images):
                self.logger.info(f"Image {i}:")
                for j, atom in enumerate(atoms):
                    self.logger.info(f"Atom {j}: {atom.symbol}, Position: {atom.position}")
        self.neb = NEB(self.images, allow_shared_calculator=True)  # NOTE: if parallelized, can't use shared calculator

        if interpolate:
            self.neb.interpolate(apply_constraint=True)
        self.set_dynamics()

    def set_dynamics(self) -> None:
        """
        Extend the set_dynamics function in the parent Job class by modifying the logger
        Redirect output to both stdout and a file
        """
        # don't call super(), neb is special
        # super().set_dynamics(name)
        self.dynamics = self.optimizer(self.neb,**self.dyn_args)
        self.dynamics.log = opt_log.__get__(self.dynamics,Optimizer)
        self.dynamics.attach(lambda : self.dyn_logger.info(self.dynamics.log()),interval=1)

    def calculate(self) -> None:
        """
        Perform the NEB Calculation
        """
        for image in self.images:  # NOTE: only works when shared calculator is enabled
            image.calc = self.potential
        max_force = self.job_params["fmax"]
        max_steps = self.job_params["max_steps"]

        steps = 0
        finished = False
        while not finished:
            self.dynamics.run(fmax=max_force,steps=1)
            steps += 1
            if steps == max_steps:
                if self.logger is not None:
                    self.logger.info('Reached NSW')
                finished = True
            # if self.get_fmax(self.curr_structure)<max_force:
            #     finished = True
        # self.create_xdatcar()
