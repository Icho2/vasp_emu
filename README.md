
# VASPyML (Name and logo TBD)

## Overview

## What we need to do:

1. [ ] Add ML potential (Open Catalyst and pyAMFF)
2. [ ] Add NEB
	-solid state NEB to optimize lattice constants
3. [ ] Add Dimer
	-solid state Dimer to optimize lattice constants
4. [ ] Add MD and INCAR-associated settings
5. [ ] Update OUTCAR (decide which output we want to emulate)
6. [ ] Make XDATCAR
7. [ ] Make OSZICAR (basically a log file of stdout)
8. [ ] Switch from ConfigParser to own implementation 
	* that will only process ones we want
	* look at TOML to see if better fit and would work without writing on
	* look to see if can just override some ConfigParser functions to make it easier
	* Make sure obey all INCAR rules (line continuation, comments, etc. . )
	* idea: just remove the check for valid options and do something else for user to signal those ignored
9. [ ] Figure out optimizers
10. [ ] Testing/examples
11. [ ] Documentation
