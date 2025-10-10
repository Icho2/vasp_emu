"""Modules"""
import os
import logging
from math import sqrt
from abc import ABC, abstractmethod

import ase
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
import ase.io
from ase.md.verlet import VelocityVerlet
from ase.optimize import BFGS, FIRE, MDMin
from ase.optimize.sciopt import SciPyFminCG
from vasp_emu.opt.sdlbfgs import SDLBFGS

class Job(ABC):
    """ 
    A class that represents a Job that can be run by the VASP Emulator
    
    Attributes:
        job_name (str) : name of the job, used for printing
        poscar (ase.Atoms) : the initial structure of the job
        job_params (dict): parameters needed for the job to run (e.g. max_steps)
        dyn_args (dict): arguments that need to be passed to the dynamics/optimizer object
        dynamics (Dynamics): the actual Dynamics/optimizer object used to guide the job
        dyn_logger (logging.Logger) : logger that is used by dynamics/optimizer object
        logger(logging.Logger) : logger that is used for job informatoin
        potential : The potential used to run the dynamics
    """
    def __init__(self, structure:ase.Atoms, dyn_name:str, dyn_args:dict, job_params:dict,
                 logger:logging.Logger=None) -> None:
        """
        Construct a Job object
        
        Arguments:
            structure (ase.Atoms): The initial structure of the job
            dyn_name (str): A string that indicates what dynamics/optimizer to use
            dyn_args (dict): The arguments needed to initializee the dynamics/optimizer
            job_params (dict): The arguments needed to run the job (e.g. max_steps, fmax)
            logger (logging.Logger): A logger object for the OUTCAR file
        """
        # Attributes to be set later
        self.job_name = ""
        self.poscar = structure
        self.potential = None
        self.dyn_logger = None
        # Other Attributes that
        self.logger = logger
        self.dyn_args = dyn_args
        self.job_params = job_params
        self.set_optimizer(dyn_name)

    @abstractmethod
    def calculate(self) -> None:
        """ An abstract method that must be defined by each instance of the Job class"""
    
    def set_dynamics(self) -> None:
        curr_structure = self.poscar
        self.dynamics = self.optimizer(curr_structure,**self.dyn_args)

    def set_optimizer(self, name:str) -> None:
        """
        Set the dynamics/optimizer object that will be used to guide the Job
        
        Arguments:
                name (str): Specifies the dynamics/optimizer object to be used
        """
        if name == "BFGS":
            self.optimizer = BFGS
        elif name == "CG":
            self.optimizer = SciPyFminCG
        elif name == "FIRE":
            self.optimizer = FIRE
        elif name == "MDMin" or name == "QuickMin":
            self.optimizer = MDMin
        elif name == "SD":
            raise NotImplementedError("SD")
        elif name == "SDLBFGS":
            self.optimizer = SDLBFGS
        elif name == "MD":
            self.optimizer = VelocityVerlet
        else:
            raise ValueError(f"Unknown dynamics type '{name}'")
        self.set_dynamics_logger()

    def set_dynamics_logger(self, name="dynamics.log") -> None:
        """Sets up the logger for the dynamics/optimizer object"""
        self.dyn_logger = logging.getLogger('opt_log')
        self.dyn_logger.setLevel(logging.INFO)
        file_handler = logging.FileHandler(name)
        file_handler.setFormatter(logging.Formatter('%(message)s'))
        self.dyn_logger.addHandler(file_handler)

    def set_potential(self, ptype:str, pname:str=None, model:str=None, infer:bool=False, device:str=None, seed:int=1234, use_cpu:bool=True) -> None:
        """
        Set the potential that will be used to find energy and forces
        
        Arguments:
                ptype (str): The type of potential to be used (e.g. UMA, PYAMFF)
                pname (str): The name of a specific model to be used when using UMA
                seed (int): Used for reproducability (default 1234)
                use_cpu (bool): whether to run the calculator on cpu (default True)
        """
        if ptype == "UMA":
            from fairchem.core import pretrained_mlip, FAIRChemCalculator
            import torch._dynamo
            torch._dynamo.config.suppress_errors = True
            from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
            from fairchem.core.units.mlip_unit import load_predict_unit
            import torch
            if infer:
                """Implementing inference setting for UMA"""
                settings = InferenceSettings(
                        tf32=True,
                        activation_checkpointing=False,
                        merge_mole=True,
                        compile=True,
                        wigner_cuda=True,
                        external_graph_gen=False,
                )
            else:
                settings = InferenceSettings(
                        tf32=False,
                        activation_checkpointing=False,
                        merge_mole=False,
                        compile=False,
                        wigner_cuda=False,
                        external_graph_gen=False,
                )
            predictor = pretrained_mlip.get_predict_unit(model, device=device, inference_settings=settings)
            self.potential = FAIRChemCalculator(predictor, task_name=pname)
                

        elif ptype == "VASP":
            ase_vasp_command = os.environ['ASE_VASP_COMMAND']
            os.system("vasp_std")
            #self.potential.read_incar("INCAR")
        elif ptype == 'PYAMFF':
            try:
                from pyamff.ase_calc import aseCalc
                from pyamff.ase_calc_fortran import aseCalcF
            except Exception as e:
                raise ModuleNotFoundError("PyAMFF potential cannot be used without PyAMFF") from e
            if os.path.isfile("mlff.pyamff"):
                self.potential = aseCalcF()
            elif os.path.isfile('pyamff.pt'):
                self.potential = aseCalc('pyamff.pt')
            else:
                raise FileNotFoundError("No PyAMFF potential file found.")
        elif ptype == "EMT":
            self.potential = EMT()
        else:
            raise ValueError(f"Unknown potential type '{ptype}' given")

    def get_fmax(self, atoms:ase.Atoms) -> float:
        """Returns fmax, as used by optimizers with NEB."""
        forces = atoms.get_forces()
        return sqrt((forces ** 2).sum(axis=1).max())

    def create_xdatcar(self, delete:bool=True) -> None:
        """
        Save the trajectory to a text file known as XDATCAR
        
        Arguments:
                delete (bool): whether to delete the intermediate trajectory file used 
                                (default True)
        """
        traj = ase.io.trajectory.Trajectory(self.dyn_args["trajectory"])
        ase.io.vasp.write_vasp_xdatcar("XDATCAR", traj)
        if delete:
            os.remove(self.dyn_args["trajectory"])
