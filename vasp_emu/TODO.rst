OUTCAR!!!!

bin/vasp_emu
============
- Determine which OCP models will work
- Add potential flags for OCP to CMD and INCAR

Emulator Class
==============

__init__
---------
- Print job type and parameters
- Should we give the user a choice in what they want to name the trajectory file?
- If a setting is overwritten on the cmd line, cmd takes precedence

read_incar()
-----------
- Handle case of ediffg like vasp

verify_dynamics()
-----------------
- Akksay doesn't like CG, ask about this
- SDLBFGS needs lglobal,lautoscale,invcurv,llineopt,fdstep?
- SD needs SDALPHA?

run()
-----
- Figure out a better way to handle NEB initial/final image, should we force the CMD?

clean()
-------
- Don't just remove, also backup
- Give option to remove pretrained models
- Make all output traj files = dynamics.traj or dyn.traj


Job Class
==========

set_dynamics()
-------------
- For SD, do we want to use QuasiNewton (ASE or Scipy) or create our own?
- QuickMin is not part of ASE or I can't find. Create our own?
- MD needs to be changed, there are multiple different algorithms we can employ

set_potential()
---------------
- Do we want to include potentials like EMT and EAM?
- Needs to give users a list of known potentials to choose from and then take that into account
- VASP: User-specify, the command (e.g. on TACC, not supposed to use mpirun)
- PyAMFF: Keep the setup or have user-specify, whether using pt or mlff
- EMT: If keeping, need to add in argumets to INCAR or cmd


OpJob Class
===========


MDJob Class 
===========

__init__
----------
- Actually read in temperature flags (T_init,T_final)

calculate()
----------
- example for Velocity Verlet, others would be different, might need to change things
- various types need to add them in

NEBJob Class
============
  
__init__
--------
- Figure out which is just for the example, and what's needed universally
- The loggers should totally be set to debug level instead of info
- What happens if someone gives a traj file with multiple images? I think this should be supported.
   - Command line option or just scan?
   - Also, should support images being given in directories as POSCARs (like VASP) 

set_dynamics()
--------------
-  the only thing that differs between all is the definition of opt_log and the parent class so, I should pass those in and just always call super
-  I could also just make the optimzer one default, nvm, that would fail with the MD

calculate()
-----------
- Multiple images so what should we output, per image?
- Same question about the CONTCARS and XDATCAR . . . 
- AH, the folders . . .

DimerJob Class
==============
- get it working, or well check if it's right

__init__
----------
- Add the dimer tags to INCAR and implement
  
opt_log()
----------
- There are two loggers at play, they have the following formatting:
   - #DIM:ROT: OPT-STEP ROT-STEP CURVATURE ROT-ANGLE ROT-FORCE
   - #MinModeTranslate: STEP      TIME          ENERGY    MAX-FORCE     STEPSIZE    CURVATURE  ROT-STEPS
   - from teh DimerControl and MinModeTranslate respectively
   - I suppressed DimerControl, should we bring back?

set_dynamics()
-------------
- Just need to figure out how I want to implement it . . .

calculate()
----------
- Is it specific to the example or for all cases, seems specific
- Are these specific controls (Dimer, displacement vector values) just for example, if so need to be universal
- trajectory should be passed from args
- don't create path.traj since we have dimer.
- dim_rlx should be replaced with self.dynamics.run(**)
- user-specify traj file?

General
========
- Add __str__ methods for all classes
- Documentation!!!
- logger.info("Reached NSW") isn't printing. . . not sure why
- All jobs doe the calc=self.potential, should I just move it to the Job class?
    -  No, NEB would have problems . . .
 - INCAR can't read comments on same line as tag, figure out what that's about and if fixable
    - Okay, let's be honest, totally fixable. Just find and strip. Just need to do it