""" Modules"""
import os
import sys
import time
import logging

import ase
import subprocess
from vasp_emu.job.job import Job
from vasp_emu.job.neb import NEBJob
from vasp_emu.job.dynamics import MDJob
from vasp_emu.job.dimer import DimerJob
from vasp_emu.job.optimization import OptJob
from vasp_emu.io.outcar import OutcarWriter
from vasp_emu.io.dimcar import DimcarWriter
from vasp_emu.utils.utils import get_sys_info
from vasp_emu.utils.config import ConfigClass
from icecream import ic

class Emulator():
    """ 
    A class used to emulate VASP
    
    Attributes:
        job (vasp_emu.job): an instance of the Job class
        oszicar (logging.Logger) : Logger that records all standard output
        user_settings (dict) : dictionary of settings passed from the command line
        config (dict) : a dictionary of settings read from the INCAR file
        params (dict) : a dictionary that maps INCAR settings to function arguments
        dynamics_name (str): the name of the dynamics object
        dyn_flags (dict): a dictionary that contains only the INCAR settings that are needed
                                by the specific dynamics object being used
    """
    def __init__(self, settings:dict=None) -> None:
        """
        Construct an Emulator object
        
        Arguments:
            settings (dict): A dictionary of settings passed from the command line
        """
        if settings.reset:
            self.clean()
            sys.exit()

        # Attributes that will be initialized in another function
        self.job = None
        self.oszicar:logging.Logger = None
        # OUTCAR will be set up in each job
        self.dyn_flags = {}
        # Settings passed from the incar
        self.config = self.read_incar(settings.INCAR)
        # This is for matching INCAR parameters with function arguments
        self.params = {**self.config,
                        "maxstep": self.config["maxmove"],
                        "trajectory": "dynamics.traj",
                        "alpha": self.config["alpha"],
                        "dt": self.config["timestep"],
                        "dtmax": self.config['ftimemax'],
                        "finc": self.config["ftimeinc"],
                        "fdec": self.config["ftimedec"],
                        "nmin": self.config["fnmin"],
                        "fa": self.config["fadec"],
                        "a": self.config["falpha"],
                        "astart": self.config["fastart"],
                        "damping": self.config['damping'],
                        "memory": self.config['ilbfgsmem'],
                        "isif": self.config['isif'],
                        "mdalgo": self.config['mdalgo'],
                        "timestep": self.config['potim'], # yes, this is confusing, jgwi
                        "max_steps" : self.config["nsw"], 
                        "temperature_K": self.config['tebeg'],
                        "tebeg": self.config['tebeg'],
                        "andersen_prob": self.config['andersen_prob'], # Andersen Thermostat tags
                        "fixcm": True,
                        "tdamp": self.config['nhc_period'], # Nose-hoover chain Thermostat tags
                        "tchain": self.config['nhc_nchains'],
                        "tloop": self.config['nhc_nrespa'],
                        "friction": self.config['langevin_gamma'], # Langevin Thermostat tags
                        "langevin_gamma": self.config['langevin_gamma'],
                        "langevin_gamma_l": self.config['langevin_gamma_l'],
                        "fmax": -1*self.config['ediffg'] if self.config['ediffg'] < 0 else 0.01,
                        "num_img": self.config['images'] if 'images' in self.config else 0,
<<<<<<< HEAD
                        "isif": self.config['isif']
=======
                        "gga": self.config['gga'],
                        "ediffg": self.config['ediffg'],
                        "ediff": self.config['ediff'],
                        "smass": self.config['smass'],
                        "ispin": self.config['ispin'],
                        "lreal": self.config['lreal'],
                        "prec": self.config['prec'],
                        "istart": self.config['istart'],
                        "iopt": self.config['iopt'],
                        "ichain": self.config['ichain'],
                        
>>>>>>> c1fe1a9ec4300e95bfc03a984eb3435212fbe025
        }
        self.dynamics_name = self.verify_dynamics(iopt=self.config["iopt"],
                                                    ibrion=self.config["ibrion"])


    def read_incar(self,filename:str='INCAR') -> dict:
        """
        Read the INCAR for job settings
        
        Arguments:
                filename (str): name of the INCAR file (default INCAR)
        Returns:
                config (dict): A dictionary with the INCAR settings
        """
        if not os.path.isfile(filename):
            print(f"Specified INCAR {''.join(os.path.abspath(filename))} does not exist",sys.stderr)
            sys.exit(2)
        # Only if we have successfully parsed the INCAR do we set up the logger
        self.setup_oszicar()
        self.oszicar.info("=======================================================")
        self.oszicar.info("Starting a VASP_EMU job at %s", time.strftime("%Y-%m-%d %X %Z"))
        self.oszicar.info("=======================================================")
        self.oszicar.info("Reading inputs from %s", filename)
        self.oszicar.info("=======================================================")
        config = ConfigClass()
        config.initialize(filename)
        return config.config

    def setup_oszicar(self) -> None:
        """Sets up OSZICAR file. OUTCAR will be set up in each job"""
        # root logger so it will take anything that's outputted
        # NOTE: Akksay says to remove this one
        self.oszicar = logging.getLogger()
        self.oszicar.setLevel(logging.INFO)
        stdout_handler = logging.StreamHandler(sys.stdout)
        stdout_handler.setFormatter(logging.Formatter('%(message)s'))
        self.oszicar.addHandler(stdout_handler)
        file_handler1 = logging.FileHandler("OSZICAR")
        file_handler1.setFormatter(logging.Formatter('%(message)s'))
        file_handler1.setLevel(logging.DEBUG)
        self.oszicar.addHandler(file_handler1)
        get_sys_info(self.oszicar)

    def verify_dynamics(self,iopt:int,ibrion:int) -> str:
        """
        Set the dynamics flags and name based on IOPT and IBRION
        
        Arguments:
                iopt (int): Specifies the optimizer to be used
                ibrion (int): Specifies optimizer if IOPT = 0, otherwise MD simulation
        Returns:
                Name of the dynamics that will be used to be passed to the Job
        """
        if iopt == 0: # we want to use the same settings as vasp
            if ibrion == 0:
                self.use_md = True
                if self.params["mdalgo"] == 1 and self.params["andersen_prob"] == 0.0: # NVE from Andersen Thermostat tags 
                    keys = ["trajectory","timestep"]
                elif self.params["mdalgo"] == 1 and self.params["isif"] == 2 and self.params["andersen_prob"] != 0.0: # Canonical NVT Ensemble with Andsersen Thermostat  
                    keys = ["trajectory","timestep","temperature_K", "andersen_prob", "fixcm"]
                elif self.params["mdalgo"] == 2 and self.params["isif"] == 2: # Canonical NVT Ensemble with Nose-Hoover Thermostat
                    keys = ["trajectory","timestep","temperature_K", "tdamp", "tchain", "tloop"]
                elif self.params["mdalgo"] == 3 and self.params["isif"] == 2: # Canonical NVT Ensemble with Langevin Thermostat
                    keys = ["trajectory", "timestep", "temperature_K", "friction"]
                elif self.params["mdalgo"] == 4 and self.params["isif"] == 2: # Canonical NVT Ensemble with Nose-Hoover Chain Thermostat  
                    keys = ["trajectory","timestep","temperature_K", "tdamp", "tchain", "tloop"]
                elif self.params["mdalgo"] == 5 and self.params["isif"] == 2: # Canonical NVT Ensemble with Canonical Sampling through Velocity Rescaling Thermostat  
                    self.logger.error("CSVR Thermostat is not yet implemented. Please try another Thermostat.")
                    sys.exit()
                elif self.params["mdalgo"] == 13 and self.params["isif"] == 2: # Canonical NVT Ensemble with Multiple Andsersen Thermostat  
                    self.logger.error("Multiple Andersen Thermostat is not yet implemented. Please try another Thermostat.")
                    sys.exit()
                self.dyn_flags = {key: self.params[key] for key in keys}
                self.dyn_flags['timestep'] *= ase.units.fs # convert to femtoseconds
                return "MD"
            if ibrion == 1:
                keys = ["trajectory","maxstep","alpha"]
                self.dyn_flags = {key: self.params[key] for key in keys}
                return "BFGS"
            if ibrion == 2:
                keys = ["trajectory"]
                self.dyn_flags = {key: self.params[key] for key in keys}
                return "CG"
            if ibrion == 3:
                keys = ["trajectory","maxstep","dt"]
                self.dyn_flags = {key: self.params[key] for key in keys}
                return "MDMin"
            raise ValueError(f"Unsupported IBRION flag: {ibrion}")
        if iopt == 1:
            if ibrion == 1:
                raise ValueError("IOPT and IBRION cannot both be 1.")
            params = ["trajectory","maxstep","damping"]
            self.dyn_flags = {key: self.params[key] for key in params}
            return "SDLBFGS"
        if iopt ==  2:
            return "CG"
        if iopt == 3:
            return "QuickMin"
        if iopt == 4:
            return "SD"
        if iopt == 7:
            params = ["trajectory","maxstep","dtmax","finc","fdec","astart","fa","a"]
            self.dyn_flags = {key: self.params[key] for key in params}
            return "FIRE"
        raise ValueError(f"Unsupported IOPT flag: {iopt}")

    def run(self) -> None:
        """Run the emulator"""
<<<<<<< HEAD
        job_params = {key: self.params[key] for key in ["max_steps","fmax", "isif"]}
=======
        job_params = self.params
>>>>>>> c1fe1a9ec4300e95bfc03a984eb3435212fbe025

        outcar_writer = OutcarWriter(isif = self.params['isif'])
        if self.config['ichain'] == 0:  # NEB
            job_params["num_img"] = self.params["num_img"]
            job_params['lclimb'] = self.config['lclimb']
            job_params['spring'] = self.config['spring']
            self.job = NEBJob(
                        structure = None,
                        dyn_name = self.dynamics_name,
                        dyn_args = self.dyn_flags,
                        job_params = job_params,
                        )
        elif self.config['ichain'] == 2:
            if (self.config['ibrion'] == 3) and (self.config['potim'] == 0):
                structure = ase.io.read("POSCAR")
                dimcar_writer = DimcarWriter()
                #outcar_writer = OutcarWriter()
                self.job = DimerJob(
                        structure = structure,
                        dyn_name = self.dynamics_name,
                        dyn_args = self.dyn_flags,
                        job_params = job_params,
                        outcar_writer = outcar_writer,
                        dimcar_writer = dimcar_writer,
                        )
            else:
                raise AttributeError("To run dimer, INCAR must include IBRION=3 and POTIM=0")
        elif self.config['ichain'] != -1:  # default value (MD or optimization)
            raise ValueError(f"Unsupported ICHAIN flag: {self.config['ichain']}\n"
                             "Only ICHAIN=0 (NEB) and ICHAIN=2 (Dimer) are currently supported in vasp_emu.")
        elif self.config["ibrion"] == 0:
            job_params["mdalgo"] = self.params["mdalgo"]
            job_params["isif"] = self.params["isif"]
            job_params["andersen_prob"] = self.params["andersen_prob"]
            job_params["tebeg"] = self.params["temperature_K"]
            structure = ase.io.read("POSCAR")
            self.job = MDJob(
                       structure = structure,
                       dyn_name = self.dynamics_name,
                       dyn_args = self.dyn_flags,
                       job_params = job_params,
                       outcar_writer = outcar_writer,
                       )
        else:
            structure = ase.io.read("POSCAR")
            #outcar_writer = OutcarWriter()
            self.job = OptJob(
                       structure = structure,
                       dyn_name = self.dynamics_name,
                       dyn_args = self.dyn_flags,
                       job_params = job_params,
                       outcar_writer = outcar_writer,
                       )
        
        if job_params['custom_model'] != 'None':
            model = os.getenv("PWD") + "/" + self.config['custom_model']
        else:
            model = self.config["umamodel"]

        self.job.set_potential(ptype=self.config["potential"], pname=self.UMA_potential(), model=model, infer=self.config["inference"], device=self.config["device"])
        self.job.calculate()
        if job_params['ml_helper'] != 'None' and job_params['initial_nsw'] != 0:
            # Here we call the second potential defined by ml_helper and then run to NSW defined by the NSW tag.
            structure = ase.io.read("CONTCAR") # CONTCAR because we are picking up where we left off.
            if job_params['ml_helper'] == 'UMA' and job_params['finetune'] == True:  # I want to finetune if we have the UMA potential
                from vasp_emu.utils.finetune_data import train_test_splits
                train_test_splits()# First make the train test splits
                subprocess.run(f"bash $HOME/vasp_emu/vasp_emu/utils/finetune.sh", shell=True, check=True)
            pname = os.getenv("PWD") + "/finetuned_model/checkpoints/final/inference_ckpt.pt"

            self.job.set_potential(ptype=self.config['ml_helper'], pname=self.config['umapot'], model=self.config["umamodel"], infer=self.config["inference"], device=self.config["device"])
            self.job.calculate() # This works. Now lets finetune before 

    def clean(self) -> None:
        """Clean up files generated by the emulator"""
        for file in ["CONTCAR","OSZICAR","OUTCAR","XDATCAR","dynamics.log","dynamics.traj"]:
            if os.path.isfile(file):
                os.remove(file)

        # remove neb files if they exist
        for i in range(1000):  # assume no one in their right mind will use more than 1000 images
            dir_name = NEBJob.neb_dir_name(i)
            if os.path.isdir(dir_name):
                for file in ["CONTCAR", "OUTCAR"]:
                    if os.path.isfile(os.path.join(dir_name, file)):
                        os.remove(os.path.join(dir_name, file))
            else:
                break
    
    def UMA_potential(self):
        if self.config['potential'] != 'UMA':
            return None
        return self.config['umapot'].lower()
