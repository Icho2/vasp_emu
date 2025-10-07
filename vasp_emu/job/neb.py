"""Modules"""
import sys, os
import time
from math import sqrt

import numpy as np

import ase
from ase.mep.neb import NEB, NEBState
from ase.optimize.optimize import Optimizer, OptimizableAtoms
from ase.geometry import find_mic
from vasp_emu.job.job import Job
from vasp_emu.io.outcar import OutcarWriter
from ase.io import write 
from ase.io.vasp import write_vasp_xdatcar

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


class VaspNEB(NEB):
    """
    A custom NEB class that intercepts intermediate values (tangents, projections)
    from ASE's internal calculations for use in the OUTCAR writer.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # This dictionary will store the stats for each image
        self.intermediate_stats = {}

    def get_forces(self):
        # --- Part 1: Replicate the data setup from ASE's BaseNEB get_forces() ---
        images = self.images
        real_forces_array = np.empty(((self.nimages - 2), self.natoms, 3))
        energies = np.empty(self.nimages)

        if self.method != 'aseneb':
            energies[0] = images[0].get_potential_energy()
            energies[-1] = images[-1].get_potential_energy()

        for i in range(1, self.nimages - 1):
            real_forces_array[i - 1] = images[i].get_forces()
            energies[i] = images[i].get_potential_energy()

        # --- Part 2: Intercept Intermediate Values ---
        state = NEBState(self, images, energies)
        self.imax = state.imax
        self.emax = state.emax # Store the barrier height

        spring1 = state.spring(0)
        for i in range(1, self.nimages - 1):
            spring2 = state.spring(i)
            
            # Use ASE's own method to get the tangent it will use
            tangent = self.neb_method.get_tangent(state, spring1, spring2, i)
            
            # Get the projection of the real force onto the tangent
            tangential_force = np.vdot(real_forces_array[i - 1], tangent)

            # Get the projection of the spring force
            proj_spring = (spring2.nt * spring2.k - spring1.nt * spring1.k)

            # Store the values we need for the OUTCAR
            self.intermediate_stats[i] = {
                'proj_real': tangential_force,
                'proj_spring': proj_spring
            }
            spring1 = spring2

        # --- Part 3: Call the original get_forces method ---
        # Now that we've stored the data we need, we let the original method
        # run to get the final effective forces for the optimizer.
        return super().get_forces()


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
            #if i != 0 and i != n_images + 1:  # don't log the first and last images
            self.outcar_writers[i] = OutcarWriter(i_dir)

        #self.neb = VaspNEB(self.images, allow_shared_calculator=True,)  # NOTE: if parallelized, can't use shared calculator
        self.neb = VaspNEB(self.images, allow_shared_calculator=True, 
                           climb=self.job_params.get('LCLIMB', True),
                           k=abs(self.job_params.get('SPRING', -5.0)),
                           method='improvedtangent') # Ensures VASP-like tangent
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


    def calculate_neb_stats_for_image(self, image_index: int) -> dict:
        """
        Retrieves intermediate NEB statistics for a given image for OUTCAR writing.
        """
        i = image_index
        stats = self.neb.intermediate_stats.get(i, {})

        # Distances and angles still need to be calculated here
        vec_prev, _ = find_mic(self.images[i-1].get_positions() - self.images[i].get_positions(), self.images[i].get_cell(), self.images[i].pbc)
        vec_next, _ = find_mic(self.images[i+1].get_positions() - self.images[i].get_positions(), self.images[i].get_cell(), self.images[i].pbc)
        
        dist_prev = np.linalg.norm(vec_prev)
        dist_next = np.linalg.norm(vec_next)

        angle = 0.0
        if dist_prev > 1e-6 and dist_next > 1e-6:
             dot = np.dot((-vec_prev).ravel(), vec_next.ravel()) / (dist_prev * dist_next)
             angle = np.arccos(np.clip(dot, -1, 1)) * 180 / np.pi

        return {
            "dist_prev": dist_prev,
            "dist_next": dist_next,
            "angle": angle,
            "proj_spring": stats.get('proj_spring', 0.0),
            "proj_real": stats.get('proj_real', 0.0)
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
        img_atoms = {}
        while not finished:
            self.dynamics.run(fmax=max_force, steps=1)
            effective_forces = self.neb.get_forces().reshape(self.neb.nimages - 2, -1, 3)
            finished = self.dynamics.converged(effective_forces.ravel())
            if step >= max_steps:
                finished = True
            step += 1

            for i, image in enumerate(self.images[1:-1], start=1):
                # Write CONTCAR
                i_dir = NEBJob.neb_dir_name(i)
                contcar = os.path.join(i_dir, "CONTCAR")
                write(contcar, image, append=False)
                
                # Save the images into img_atoms dict
                try:
                    img_atoms[i_dir].append(image)
                except:
                    img_atoms[i_dir] = [image]

                # 1. Get all common stats (from base Job class)
                image_forces = effective_forces[i-1]
                step_data = self.get_step_data(image, image_forces, calc_stress = True)
                
                # 2. Get all NEB-specific stats (from this class)
                neb_stats_dict = self.calculate_neb_stats_for_image(i)
                
                # 3. Pass ALL data to the new write_step method
                self.outcar_writers[i].write_step(
                    step, 
                    step_data=step_data, 
                    neb_stats=neb_stats_dict
                )

        # TODO: figure out how to create XDATCAR for each image
            for i_dir in img_atoms:
                xdatcar = os.path.join(i_dir, "XDATCAR")
                write_vasp_xdatcar(xdatcar, img_atoms[i_dir])



    @staticmethod
    def neb_dir_name(i:int) -> str:
        """
        Convert image number to directory name for NEB calculations
        (e.g. 0 -> 00, 1 -> 01, ..., 10 -> 10)
        """
        return str(i) if i > 9 else "0" + str(i)
