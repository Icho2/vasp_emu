"""Modules"""
import time
from math import sqrt

import ase
from ase import units
from ase.optimize.optimize import Dynamics
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from vasp_emu.job.job import Job
from icecream import ic

def opt_log(self, forces=None) -> str:
    """
        Redefine the behavior of the log function for Dynamics. 
        
        Arguments:
            forces (list of floats) : The forces that correspond to the optimizable object
        Returns:
            Message to be sent to the logger
    """
    if forces is None:
        forces = self.optimizable.atoms.get_forces()
    if forces.ndim == 1:
        forces = forces.reshape(-1,3)
    fmax = sqrt((forces ** 2).sum(axis=1).max())
    #epot = self.optimizable.atoms.get_potential_energy() / len(self.optimizable.atoms)
    epot = self.optimizable.atoms.get_potential_energy()
    #ekin = self.atoms.get_kinetic_energy() / len(self.optimizable.atoms)
    ekin = self.atoms.get_kinetic_energy() / len(self.optimizable.atoms)
    etot = epot+ekin
#    eext = self.dynamics.get_extended_potential_energy()
    t = time.localtime()
    name = self.__class__.__name__
    # everything above this line exactly matches the Optimizer.log()
    msg = ""
    if self.nsteps == 0:
        msg += "All energies are eV per atoms\n"
        args = (" " * len(name), "Step", "Time", " Epot ", " Ekin ",  " Etot ", "fmax", "Temp(K)")
        msg += "%s  %4s %10s %14s %12s %12s %10s %12s\n" % args

    args = (name, self.nsteps, t[3], t[4], t[5], epot, ekin, etot, fmax, ekin / (1.5 * units.kB))
    msg += "%s:  %3d     %02d:%02d:%02d %12.6f %12.6f %12.6f %12.6f %10.2f" % args
    return msg

class MDJob(Job):
    """ 
    An instance of the Job class used to run molecular dynamics
    
    Attributes:
        job_name (str): name of the job
    """ 
    def __init__(self,T_init:int=300,T_final:int=None,**kwargs):
        """
        Construct an OptJob
        
        Arguments:
            init_temp (int) : initial temperature of the MD simulation
            final_temp (int) : final temperature of the MD simulation
            Please refer to the parent Job Class
        """
        super().__init__(**kwargs) # always first
        self.job_name = "molecular-dynamics"
        self.init_temp = T_init
        self.final_temp = T_final if T_final is not None else T_init
        self.set_dynamics() # always last

    def set_dynamics(self) -> None:
        """
        Extend the set_dynamics function in the parent Job class by modifying the logger
        Redirect output to both stdout and a file
        """
        super().set_dynamics()
        self.dynamics.log = opt_log.__get__(self.dynamics,Dynamics)
        self.dynamics.attach(lambda : self.dyn_logger.info(self.dynamics.log()),interval=1)

    def calculate(self):
        """
        Perform the MD simulation
        """
        curr_structure = self.poscar
        max_steps = self.job_params["max_steps"]
        curr_structure.calc = self.potential
        steps = 0
        finished = False
        MaxwellBoltzmannDistribution(curr_structure,temperature_K=self.job_params["tebeg"])
        '''Here I will add my sauce'''
        while not finished:
            self.dynamics.run(steps=1)
            ic(self.dynamics.get_conserved_energy())
            ic(self.dynamics._thermostat.get_thermostat_energy())
            self.logger.info(f'U: {curr_structure.get_potential_energy()}   ' + \
                                f'fmax: {self.get_fmax(curr_structure)}')
            # CONTCAR should be written after each step, used to restart jobs
            ase.io.write('CONTCAR',curr_structure,format='vasp')
            steps += 1
            if steps == self.job_params['initial_nsw'] and self.job_params['ml_helper'] != 'None': # Stop MD and Fintune your model or make it ! 
                # For now I am not finetuning because I need to figure out how to switch potentials first
                # Will not switch potentials here, instead I will end the calculation and start another with a new potential in emulate.py
                self.logger.info('Finished initial run, now running the second potential')
                finished = True
            if steps == max_steps:
                self.logger.info('Reached NSW')
                finished = True
        self.create_xdatcar(False)
