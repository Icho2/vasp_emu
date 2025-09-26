"""Modules"""
import sys, os
import time
from math import sqrt

import numpy as np

import ase
from ase.mep import NEB
from ase.optimize.optimize import Optimizer, OptimizableAtoms
from ase.geometry import find_mic
from vasp_emu.job.job import Job
from vasp_emu.io.outcar import OutcarWriter
from ase.io import write

# In vasp_emu/job/neb.py

def opt_log(self, forces=None) -> str:
    """
    Redefine the behavior of the log function for Optimizers when used with NEB.
    """
    # For NEB, self.optimizable is a NEBOptimizable object.
    # We need to get the forces and energy from the NEB object it contains.
    neb = self.optimizable.neb
    
    if forces is None:
        forces = neb.get_forces() # Get effective forces from the NEB object

    # The fmax is the maximum effective force across all images
    fmax = sqrt((forces.reshape(-1, 3) ** 2).sum(axis=1).max())
    
    # For NEB, get_potential_energy() returns the highest energy of the band (emax)
    e = neb.get_potential_energy()
    
    t = time.localtime()
    name = self.__class__.__name__
    
    msg = ""
    if self.nsteps == 0:
        msg += "=======================================================\n"
        # The energy here is the barrier height, not a total energy.
        msg += f'{" " * len(name)}  {"Step":4s} {"Time":>9s} {"Barrier":>13s}  {"fmax":>10s}\n'

    msg += f"{name}:  {self.nsteps:3d}    {t[3]:02d}:{t[4]:02d}:{t[5]:02d} {e:12.6f} {fmax:12.6f}"
    return msg


class NEBJob(Job):
    """ 
    An instance of the Job class used to perform nudged elastic band method
    
    Attributes:
        job_name (str): name of the job
        images (list of Atoms): list of images in the NEB calculation
                                     self.poscar is not used
        outcar_writers (list of OutcarWriter): list of loggers for each image
                                               self.outcar_writer is not used
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
        self.outcar_writers = dict()  # dict of outcar loggers for each image
        for i in range(n_images+2):
            i_dir = NEBJob.neb_dir_name(i)
            poscar = os.path.join(i_dir, "POSCAR")
            curr_structure = ase.io.read(poscar)
            self.images.append(curr_structure)
            if i != 0 and i != n_images + 1:  # don't log the first and last images
                self.outcar_writers[i] = OutcarWriter(i_dir)

        self.neb = NEB(self.images, allow_shared_calculator=True)  # NOTE: if parallelized, can't use shared calculator
        self.set_dynamics()
        print(self.dyn_args["trajectory"])

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


    def calculate_neb_stats_for_image(self, image_index: int) -> dict:
        """
        Calculates all VASP-style NEB stats for a given image.
        """
        i = image_index
        
        # Get the relevant images from the list
        image_current = self.images[i]
        image_prev = self.images[i - 1]
        image_next = self.images[i + 1]

        # --- 1. CORRECTLY Calculate displacement vectors (dR) and distances ---
        
        # Get vector difference for previous image, accounting for periodic boundaries
        pos_diff_prev = image_prev.get_positions() - image_current.get_positions()
        dR_prev, _ = find_mic(pos_diff_prev, image_current.get_cell(), image_current.pbc)
        dist_prev = np.linalg.norm(dR_prev)

        # Get vector difference for next image, accounting for periodic boundaries
        pos_diff_next = image_next.get_positions() - image_current.get_positions()
        dR_next, _ = find_mic(pos_diff_next, image_current.get_cell(), image_current.pbc)
        dist_next = np.linalg.norm(dR_next)
        
        # --- 2. Determine the tangent (replicating VASP logic) ---
        # This logic chooses the tangent based on the energy landscape
        E_current = image_current.get_potential_energy()
        E_prev = image_prev.get_potential_energy()
        E_next = image_next.get_potential_energy()

        # Note: VASP uses -dR_prev as the forward vector from prev->current
        tangent_prev = -dR_prev 
        tangent_next = dR_next

        if E_prev > E_current < E_next: # At a minimum
            tangent = tangent_next * abs(E_next - E_current) + tangent_prev * abs(E_current - E_prev)
        elif E_prev < E_current > E_next: # At a maximum (or saddle)
            tangent = tangent_next * abs(E_next - E_current) + tangent_prev * abs(E_current - E_prev)
        elif E_next > E_prev:
            tangent = tangent_next
        else:
            tangent = tangent_prev
        
        # Normalize the tangent
        tangent_norm = np.linalg.norm(tangent)
        if tangent_norm > 0:
            tangent /= tangent_norm

        # --- 3. Calculate angle between vectors ---
        if dist_prev > 0 and dist_next > 0:
            unit_dR_prev = -dR_prev / dist_prev # Pointing forward
            unit_dR_next = dR_next / dist_next  # Pointing forward
            dot_product = np.clip(np.sum(unit_dR_prev * unit_dR_next), -1.0, 1.0)
            angle = np.arccos(dot_product) * 180.0 / np.pi
        else:
            angle = 0.0

        # --- 4. Calculate projections ---
        spring_const = abs(self.job_params.get("SPRING", -5.0))
        proj_spring = spring_const * (dist_next - dist_prev)
        
        force_real = image_current.get_forces() 
        proj_real = np.sum(force_real * tangent)
        
        return {
            "dist_prev": dist_prev,
            "dist_next": dist_next,
            "angle": angle,
            "proj_spring": proj_spring,
            "proj_real": proj_real
        }


    def get_neb_header_params(self) -> dict:
        """Helper to extract NEB header params from job_params."""
        # TODO: Pull these from your INCAR parser (which should populate self.job_params)
        return {
            "spring": self.job_params.get("SPRING", -5.0),
            "lclimb": self.job_params.get("LCLIMB", True),
            "efirst": self.images[0].get_potential_energy(), # Requires calc on first image
            "elast": self.images[-1].get_potential_energy() # Requires calc on last image
        }

    def calculate(self) -> None:
        """
        Perform the NEB Calculation (now with new logging)
        """
        for image in self.images:
            image.calc = self.potential
            
        # Ensure endpoints have energies for the header
        #if self.images[0].calc: self.images[0].get_potential_energy()
        #if self.images[-1].calc: self.images[-1].get_potential_energy()
            
        max_force = self.job_params["fmax"]
        max_steps = self.job_params["max_steps"]

        # write OUTCAR headers
        neb_header_params = self.get_neb_header_params()
        for i, image in enumerate(self.images[1:-1], start=1):  # skip first and last images
            # Pass the NEB-specific header params
            self.outcar_writers[i].write_header(image, neb_params=neb_header_params)

        step = 1
        finished = False
        #import copy
        #previous_images = [copy.deepcopy(img) for img in self.images]
        while not finished:
            self.dynamics.run(fmax=max_force, steps=1)
            finished = self.dynamics.converged(self.neb.get_forces().ravel())

            for i, image in enumerate(self.images[1:-1], start=1):
                # Write CONTCAR
                i_dir = NEBJob.neb_dir_name(i)
                contcar = os.path.join(i_dir, "CONTCAR")
                write(contcar, image, append=False)
                
                # 1. Get all common stats (from base Job class)
                step_data = self.get_step_data(image)
                
                # 2. Get all NEB-specific stats (from this class)
                neb_stats_dict = self.calculate_neb_stats_for_image(i)
                
                # 3. Pass ALL data to the new write_step method
                self.outcar_writers[i].write_step(
                    step, 
                    step_data=step_data, 
                    neb_stats=neb_stats_dict
                )
                # --- END MODIFIED SECTION ---

            if step >= max_steps:
                finished = True
            step += 1
            #previous_images = [copy.deepcopy(img) for img in self.images]

        # TODO: figure out how to create XDATCAR for each image


    @staticmethod
    def neb_dir_name(i:int) -> str:
        """
        Convert image number to directory name for NEB calculations
        (e.g. 0 -> 00, 1 -> 01, ..., 10 -> 10)
        """
        return str(i) if i > 9 else "0" + str(i)