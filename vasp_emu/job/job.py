"""Modules"""
import os
import logging
from math import sqrt
from abc import ABC, abstractmethod

import numpy as np
import ase
import ase.md
from ase.md.nose_hoover_chain import NoseHooverChainNVT 
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase.md.andersen import Andersen  
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp
import ase.io
from ase.optimize import BFGS, FIRE, MDMin
from ase.optimize.sciopt import SciPyFminCG
from vasp_emu.opt.sdlbfgs import SDLBFGS
from vasp_emu.io import outcar, dimcar

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
        outcar_writer (outcar.OutcarWriter) : logger that is used for job informatoin
        potential : The potential used to run the dynamics
    """
    def __init__(self, structure:ase.Atoms, dyn_name:str, dyn_args:dict, job_params:dict,
                 outcar_writer:outcar.OutcarWriter=None,
                 dimcar_writer:dimcar.DimcarWriter=None) -> None:
        """
        Construct a Job object
        
        Arguments:
            structure (ase.Atoms): The initial structure of the job
            dyn_name (str): A string that indicates what dynamics/optimizer to use
            dyn_args (dict): The arguments needed to initializee the dynamics/optimizer
            job_params (dict): The arguments needed to run the job (e.g. max_steps, fmax)
            outcar_writer (outcar.OutcarWriter): A logger object for the OUTCAR file
        """
        # Attributes to be set later
        self.job_name = ""
        self.poscar = structure
        self.potential = None
        self.dyn_logger = None
        # Other Attributes that
        self.outcar_writer = outcar_writer
        self.dimcar_writer = dimcar_writer
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
            if self.job_params["mdalgo"] == 1 and self.job_params["isif"] == 2 and self.job_params["andersen_prob"] != 0.0: # Canonical NVT Ensemble with Andsersen Thermostat  
                self.optimizer = Andersen
            elif self.job_params["mdalgo"] == 2 and self.job_params["isif"] == 2: # Canonical NVT Ensemble with Nose-Hoover Thermostat  
                self.optimizer = NoseHooverChainNVT # Is the chain version the same as the standard?
            elif self.job_params["mdalgo"] == 3 and self.job_params["isif"] == 2: # Canonical NVT Ensemble with Langevin Thermostat  
                self.optimizer = Langevin
            elif self.job_params["mdalgo"] == 4 and self.job_params["isif"] == 2: # Canonical NVT Ensemble with Nose-Hoover Chain Thermostat  
                self.optimizer = NoseHooverChainNVT
            elif self.job_params["mdalgo"] == 5 and self.job_params["isif"] == 2: # Canonical NVT Ensemble with Andsersen Thermostat  
                self.logger.error("CSVR Thermostat is not yet implemented. Please try another Thermostat.")
                self.optimizer = None
            elif self.job_params["mdalgo"] == 13 and self.job_params["isif"] == 2: # Canonical NVT Ensemble with Multiple Andsersen Thermostat  
                self.logger.error("Multiple Andersen Thermostat is not yet implemented. Please try another Thermostat.")
                self.optimizer = None
            elif self.job_params["mdalgo"] == 1 and self.job_params["andersen_prob"] == 0.0:
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
            model_name = {"S":"uma-s-1p1", "M":"uma-m-1p1"}
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
            if self.job_params['custom_model'] != 'None':
                predictor = pretrained_mlip.load_predict_unit(model, device=device, inference_settings=settings)
            else:
                predictor = pretrained_mlip.get_predict_unit(model_name[model], device=device, inference_settings=settings)
            self.potential = FAIRChemCalculator(predictor, task_name=pname)
                

        elif ptype == 'VASP':            
            # Make a copy of a original INCAR because this VASP command overwrites it.
            import subprocess
            executable = os.environ["ASE_VASP_COMMAND"] 
            self.potential= Vasp(command=executable,
                                    restart=False, directory='vasp_run',
                                     xc = self.job_params['gga'],
                                     nsw = self.job_params['nsw'],
                                     ediffg = self.job_params['ediffg'],
                                     ediff = self.job_params['ediff'],
                                     encut = self.job_params['encut'],
                                     smass = self.job_params['smass'],
                                     langevin_gamma = self.job_params['langevin_gamma'],
                                     tebeg = self.job_params['tebeg'],
                                     langevin_gamma_l = self.job_params['langevin_gamma_l'],
                                     ispin = self.job_params['ispin'],
                                     lreal = self.job_params['lreal'],
                                     andersen_prob = self.job_params['andersen_prob'],
                                     prec = self.job_params['prec'],
                                     istart = self.job_params['istart'],
                                     isif = self.job_params['isif'],
                                     iopt = self.job_params['iopt'],
                                     ichain = self.job_params['ichain'],
                                     images = self.job_params['images'],
                                     ibrion = self.job_params['ibrion'],
                                     damping = self.job_params['damping'],
                                     nblock  = self.job_params['nblock']
                                     )

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


    def create_xdatcar(self, delete:bool=False) -> None:
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


    def get_stress(self, atoms: ase.Atoms) -> tuple:
        """
        Returns the stress tensor in voight form, stress tensor trace, and trace / sqrt(3)
        all in a tuple.
        ISIF == 0: return None
        ISIF == 0: return trace and trace / sqrt(3) with tensor =  None
        ISIF > 0: return all 3

        Args:
            atoms (ase.Atoms): the Atoms object for any given structure
        """
        if self.job_params["isif"] == 0:
            return (None, None, None)
        elif self.job_params["isif"] > 0:
            try:
                stress = atoms.get_stress()
            except:
                # TODO above gets you stress if the calculator calculates the stress
                # otherwise we have to do it manually below
                # right now return a value of -100
                stress = [-100 for i in range(6)]
        else:
            raise ValueError('The ISIF tag can only be positive integers up to 8')
        
        trace = np.sqrt(stress[0] ** 2 + stress[1] ** 2 + stress[2] ** 2)
        dim = trace / np.sqrt(3)

        if self.job_params["isif"] == 1:
            return (None, trace, dim)
        elif self.job_params["isif"] > 1:
            return(stress, trace, dim)
            
    def get_step_data(self, atoms: ase.Atoms, forces: np.ndarray) -> dict:
        """
        Gathers ALL data needed for an OUTCAR step into a single dictionary.

        Args:
            atoms (ase.Atoms): The atoms object for the current step.
            forces (np.ndarray): The forces to be used for statistics.
                                 This can be TRUE forces or EFFECTIVE forces
                                 for NEBJob.
                                 Dimensions (n_atoms, 3).
        """
        # Basic data
        positions = atoms.get_positions()
        energy = atoms.get_potential_energy()
        volume = atoms.get_volume()

        # Derived stats calculated from the provided 'forces' array
        fmax_atom = np.sqrt((forces ** 2).sum(axis=1).max())
        f_rms = np.sqrt(np.mean(forces**2))

        # Placeholder/Parameter-based stats
        stress_tensor, stress_total, stress_dim = self.get_stress(atoms)
        n_electrons = self.job_params.get('NELECT', int(np.sum(atoms.get_atomic_numbers()))) # either get NELECT
        # or sum the atomic numbers to get the total electrons
        magnetization = self.job_params.get('NUPDOWN', 0)

        return {
            "positions": positions,
            "forces": forces,  # This now correctly stores either true or effective
            "energy": energy,
            "fmax_atom": fmax_atom,
            "f_rms": f_rms,
            "stress_total": stress_total,
            "stress_dim": stress_dim,
            "stress_tensor": stress_tensor, # voight form
            "volume": volume,
            "n_electrons": n_electrons,
            "magnetization": magnetization
        }


#    def get_fmax(self, atoms:ase.Atoms) -> float:
#        """Returns fmax, as used by optimizers."""
#        forces = atoms.get_forces()
#        return sqrt((forces ** 2).sum(axis=1).max())
