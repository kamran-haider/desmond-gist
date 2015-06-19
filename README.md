## desmond-gist
**Implementation of grid inhomogeneous solvation theory in python.**
GIST provides spatially resolved contribution to free energy of solvation (with breakdowon of entropy and enthalpy) in systems of biomolecular and pharmaceutical interest. This version of GIST contains additional features that provide information on local structue of water in protein active sites, which can be valuable for ligand design. This program requires a working Desmond/Maestro installation and is intended mainly for Desmond generated trajectories. 

##Installation
GIST is implemented as a python script with an associated C module (_gistcalc.c) that needs to be compiled. The compilation instructions in the C module can
 the main script _run_gist.py_ can be run using $SCHRODINGER/run
