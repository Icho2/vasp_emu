
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
	* Look at TOML to see if it's a better fit and would work without writing on
	* Look to see if we can override some ConfigParser functions to make it easier
	* Make sure to obey all INCAR rules (line continuation, comments, semicolons indicating secondary line, etc. . )
	* idea: remove the check for valid options and do something else for the user to signal those ignored
9. [ ] Figure out optimizers
10. [ ] Testing/examples
11. [ ] Documentation
