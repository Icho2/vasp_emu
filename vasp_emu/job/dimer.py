"""Modules"""
import time
from math import sqrt
import numpy as np

import ase
from ase.constraints import FixAtoms
from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate  
from ase.mep.dimer import normalize, rotate_vectors, DimerEigenmodeSearch

from vasp_emu.job.job import Job


def opt_log(self, forces=None):
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

class DimerJob(Job):
    """ 
        An instance of the Job class used to run dimer
        
        Attributes:
            job_name (str): name of the job
    """
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.job_name = "dimer"
        self.set_dynamics() # always last

    def get_atc(self, d_atoms: ase.mep.dimer.MinModeAtoms) -> tuple:
        """
        Gets the rotational angles and corresponding torques and curvatures for each step of dimer
        
        Args:
            d_atoms (ase.mep.dimer.MinModeAtoms): dimer atoms object from calculate

        Returns:
            a tuple(list) where first list are the angles, second are torques, and thrid are curvatures
        """
        dR = d_atoms.control.get_parameter('dimer_separation')
        n = d_atoms.get_eigenmode()
        pos1 = d_atoms.get_positions() + n * dR
        pos2 = d_atoms.get_positions() - n * dR
        search = DimerEigenmodeSearch(d_atoms, d_atoms.control, d_atoms.eigenmodes[0], d_atoms.eigenmodes[:0])

        # based on ase.mep.dimer.Dimereigenmode.converge_eigenmode
        stoprot = False

        # Load the relevant parameters from control
        f_rot_min = d_atoms.control.get_parameter('f_rot_min')
        f_rot_max = d_atoms.control.get_parameter('f_rot_max')
        trial_angle = d_atoms.control.get_parameter('trial_angle')
        max_num_rot = d_atoms.control.get_parameter('max_num_rot')
        extrapolate = d_atoms.control.get_parameter('extrapolate_forces')

        angles = [] 
        torques = []
        curvatures = []

        while not stoprot:
            f1 = d_atoms.get_forces(real = True, pos = pos1)
            if d_atoms.control.get_parameter('use_central_forces'):
                f2 = 2 * d_atoms.get_forces(real = True) - f1
            else:
                f2 = d_atoms.get_forces(real = True, pos = pos2)

            search.forces1 = f1
            search.forces2 = f2
            rA = normalize(search.get_rotational_force())
            rB = rotate_vectors(n, rA, trial_angle)[1]


            # calculate torquq as norm of rotational force
            torque = np.linalg.norm(search.get_rotational_force())
            curvature = d_atoms.get_curvature()
            
            # Pre rotation stop criteria
            if torque <= f_rot_min:
                stoprot = True
            else:

                c0 = np.vdot((f2 - f1), rA) / dR
                c1 = np.vdot((f2 - f1), rB) / dR

                # Calculate the Fourier coefficients
                a1 = c0 * np.cos(2 * trial_angle) - c1 / (2 * np.sin(2 * trial_angle))
                b1 = 0.5 * c0
                a0 = 2 * (c0 - a1)

                # Estimate the rotational angle
                rotangle = np.atan(b1 / a1) / 2.0

                # Make sure that you didn't find a maximum
                cmin = a0 / 2.0 + a1 * np.cos(2 * rotangle) + b1 * np.sin(2 * rotangle)
                if c0 < cmin:
                    rotangle += np.pi / 2.0
                

                # append rotangle to angles and torque to torques
                angles.append(rotangle)
                torques.append(torque)
                curvatures.append(curvature)

            # Post rotation stop criteria
            if not stoprot:
                if d_atoms.control.get_counter('rotcount') >= max_num_rot:
                    stoprot = True
                elif norm(search.get_rotational_force()) <= f_rot_max:
                    stoprot = True
        return (angles, torques, curvatures)

    def calculate(self) -> None:
        """
        Perform the dimer calculation
        """
        curr_structure = self.poscar
        mask = [atom.tag > 0 for atom in curr_structure]
        curr_structure.set_constraint(FixAtoms(mask=mask))
        curr_structure.calc = self.potential
        curr_structure.get_potential_energy()
        max_force = self.job_params["fmax"]
        max_steps = self.job_params["max_steps"]

        steps = 0
        finished = False

        self.dimcar_writer.write_dimcar_header()
        with DimerControl(initial_eigenmode_method='gauss', displacement_method='gauss', logfile=None, mask=[0,0,0,0,1]) as d_control:
            d_atoms = MinModeAtoms(curr_structure,d_control)

            # Displace the atoms
            #displacement_vector = [[0.0] * 3] * 5
            #displacement_vector[-1][1] = -0.1
            #d_atoms.displace(displacement_vector=displacement_vector)

            with MinModeTranslate(d_atoms, trajectory='dynamics.traj') as dim_rlx:
                while not finished:
                    dim_rlx.run(fmax=max_force,steps=1)
                    forces = curr_structure.get_forces()
                    fmax = sqrt((forces ** 2).sum(axis=1).max())

                    # Write step data
                    step_data = self.get_step_data(curr_structure, forces)
                    self.outcar_writer.write_step(
                        steps, 
                        step_data=step_data, 
                    )

                    # CONTCAR should be written after each step, used to restart jobs
                    ase.io.write('CONTCAR',curr_structure,format='vasp')
                    # write the CENTCAR after each step
                    ase.io.write('CENTCAR',curr_structure,format='vasp')
                    steps+=1

                    angles, torques, curvatures = self.get_atc(d_atoms) 
                    self.dimcar_writer.write_dimcar(
                            {"step": steps,
                             "fmax": fmax,
                             "torques": torques,
                             "energy": d_atoms.get_potential_energy(),
                             "curvatures": curvatures,
                             "angles": angles
                             })

                    if fmax < max_force:
                        finished = True
                        self.outcar_writer.info('Reached accuracy')
                    elif steps == max_steps:
                        self.outcar_writer.info('Reached NSW')
                        finished  = True
        
        # Write the MODECAR
        with open('MODECAR', 'w') as f1:
            atoms1 = ase.io.read('dynamics.traj', index = 0)
            atoms2 = ase.io.read('dynamics.traj', index = 1)
            diff = atoms2.get_positions() - atoms1.get_positions()
            for d in diff:
                f1.write(f"{d[0]:20.10E}{d[1]:20.10E}{d[2]:20.10E}\n")

        # self.create_xdatcar()
