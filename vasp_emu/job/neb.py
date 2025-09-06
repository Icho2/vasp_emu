"""Modules"""
import sys, os
import time
from math import sqrt

import ase
from ase.mep import NEB
from ase.optimize.optimize import Optimizer, OptimizableAtoms
from vasp_emu.job.job import Job
from vasp_emu.io.outcar import OutcarWriter
from ase.io import write

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
        images (list of Atoms): list of images in the NEB calculation
                                     self.poscar is not used
        loggers (list of OutcarWriter): list of loggers for each image
                                        self.logger is not used
        neb (NEB) : NEB object to drive the dynamics
        Please refer to the parent Job Class
    """
    def __init__(self,**kwargs):
        """
        Construct an NEBJob
        
        Arguments:
            Please refer to the parent Job Class
        """
        super().__init__(**kwargs)
        self.job_name = "nudged-elastic-band"

        # populate images (assume nebmake was already run by the user)
        self.images = []
        n_images = self.job_params['num_img']
        # the directory n_images+2 shouldn't exist
        if os.path.exists(NEBJob.neb_dir_name(n_images+2)):
            raise ValueError(f"Directory {NEBJob.neb_dir_name(n_images+2)} exists.\n"
                                "Double check the IMAGES tag in the INCAR file.\n")
        self.loggers = dict()  # dict of outcar loggers for each image
        for i in range(n_images+2):
            i_dir = NEBJob.neb_dir_name(i)
            poscar = os.path.join(i_dir, "POSCAR")
            curr_structure = ase.io.read(poscar)
            self.images.append(curr_structure)
            if i != 0 and i != n_images + 1:  # don't log the first and last images
                self.loggers[i] = OutcarWriter(i_dir)

        self.neb = NEB(self.images, allow_shared_calculator=True)  # NOTE: if parallelized, can't use shared calculator
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

        # write OUTCAR headers
        for i, image in enumerate(self.images[1:-1], start=1):  # skip first and last images
            self.loggers[i].write_header(image)

        step = 1  # VASP ionic steps start at 1
        finished = False
        while not finished:
            self.dynamics.run(fmax=max_force,steps=1)
            finished = self.dynamics.converged()

            # write CONTCAR to image directories
            for i, image in enumerate(self.images[1:-1], start=1):  # skip first and last images
                i_dir = NEBJob.neb_dir_name(i)
                contcar = os.path.join(i_dir, "CONTCAR")
                write(contcar, image, append=False)
                # write OUTCAR
                self.loggers[i].write_step(image, step)

            if step == max_steps:
                finished = True
            # if self.get_fmax(self.curr_structure)<max_force:
            #     finished = True
            step += 1
        # self.create_xdatcar()

    @staticmethod
    def neb_dir_name(i:int) -> str:
        """
        Convert image number to directory name for NEB calculations
        (e.g. 0 -> 00, 1 -> 01, ..., 10 -> 10)
        """
        return str(i) if i > 9 else "0" + str(i)