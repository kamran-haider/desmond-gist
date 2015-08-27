## desmond-gist
**Implementation of grid inhomogeneous solvation theory in python.**
GIST provides spatially resolved contribution to free energy of solvation (with breakdowon of entropy and enthalpy) in systems of biomolecular and pharmaceutical interest. This version of GIST contains additional features that provide information on local structue of water in protein active sites, which can be valuable for ligand design. This program requires a working Desmond/Maestro installation and is intended mainly for Desmond generated trajectories. 

##Installation
GIST is implemented as a python script with an associated C module (_\_gistcalc.c_) that needs to be compiled. For a desmond installation of python 2.5 (change path up to the location of desmond directory on your local machine, rest should be the same):
  
gcc -O3 -lm -fPIC -shared -I path_to_desmond/desmond/mmshare-v24012/lib/Linux-x86_64/include/python2.7 -I path_to_desmond/desmond/mmshare-v24012/lib/Linux-x86_64/lib/python2.7/site-packages/numpy/core/include/ -o _gistcalcs.so _gistcalcs.c

The compilation instructions are provided in the C module 
The main script _run\_gist.py_ should be run using $SCHRODINGER/run as it requires Schrodinger Python API. Type _$SCHRODINGER/run run\_gist.py_ --help_ for on the command line options required to run a gist calculation.
##Usage
For a general introduction to running GIST calculations, please visit: http://ambermd.org/tutorials/advanced/tutorial25/
Further details can be found in gist_tutorial.pdf. 
