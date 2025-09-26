""" Modules"""
import os
import sys
import time
import logging

import ase
from vasp_emu.job.job import Job
from vasp_emu.job.neb import NEBJob
from vasp_emu.job.dynamics import MDJob
from vasp_emu.job.dimer import DimerJob
from vasp_emu.job.optimization import OptJob
from vasp_emu.io.outcar import OutcarWriter
from vasp_emu.utils.utils import get_sys_info
from vasp_emu.utils.config import ConfigClass

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
        self.oszicar:logging.Logger = None # Akksay vetoed
        # OUTCAR will be set up in each job
        self.dyn_flags = {}
        # Settings passed from the incar
        self.config = self.read_incar(settings.INCAR)

        # This is for matching INCAR parameters with function arguments
        self.params = {
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
                        "timestep": self.config['potim'], # yes, this is confusing, jgwi
                        "max_steps" : self.config["nsw"], 
                        "fmax": -1*self.config['ediffg'] if self.config['ediffg'] < 0 else 0.01,
                        "num_img": self.config['images'] if 'images' in self.config else 0
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
                keys = ["trajectory","timestep"]
                self.dyn_flags = {key: self.params[key] for key in keys}
                self.dyn_flags["timestep"] *= ase.units.fs # convert to femtoseconds
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
        job_params = {key: self.params[key] for key in ["max_steps","fmax"]}

        if self.config['ichain'] == 0:  # NEB
            job_params["num_img"] = self.params["num_img"]
            self.job = NEBJob(
                        structure = None,
                        dyn_name = self.dynamics_name,
                        dyn_args = self.dyn_flags,
                        job_params = job_params,
                        )
        elif self.config['ichain'] == 2:
            if (self.config['ibrion'] == 3) and (self.config['potim'] == 0):
                structure = ase.io.read("POSCAR")
                outcar_writer = OutcarWriter()
                self.job = DimerJob(
                        structure = structure,
                        dyn_name = self.dynamics_name,
                        dyn_args = self.dyn_flags,
                        job_params = job_params,
                        outcar_writer = outcar_writer,
                        )
            else:
                raise AttributeError("To run dimer, INCAR must include IBRION=3 and POTIM=0")
        elif self.config['ichain'] != -1:  # default value (MD or optimization)
            raise ValueError(f"Unsupported ICHAIN flag: {self.config['ichain']}\n"
                             "Only ICHAIN=0 (NEB) and ICHAIN=2 (Dimer) are currently supported in vasp_emu.")
        elif self.config["ibrion"] == 0:
            structure = ase.io.read("POSCAR")
            outcar_writer = OutcarWriter()
            self.job = MDJob(
                       structure = structure,
                       dyn_name = self.dynamics_name,
                       dyn_args = self.dyn_flags,
                       job_params = job_params,
                       outcar_writer = outcar_writer,
                       )
        else:
            structure = ase.io.read("POSCAR")
            outcar_writer = OutcarWriter()
            self.job = OptJob(
                       structure = structure,
                       dyn_name = self.dynamics_name,
                       dyn_args = self.dyn_flags,
                       job_params = job_params,
                       outcar_writer = outcar_writer,
                       )

        self.job.set_potential(ptype=self.config["potential"], pname=self.UMA_potential(), model=self.config["umamodel"], infer=self.config["inference"], device=self.config["device"])
        self.job.calculate()

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
