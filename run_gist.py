# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:03:09 2015

@author: kamran
"""


#*********************************************************************************************#
import time
import numpy as np
from schrodinger import structure
from optparse import OptionParser
from gist_desmond_composite import Gist

parser = OptionParser()
parser.add_option("-i", "--input_cms", dest="cmsname", type="string", help="Input CMS file")
parser.add_option("-t", "--input_trajectory", dest="trjname", type="string", help="Input trajectory directory")
parser.add_option("-l", "--clust_center", dest="clust_center_file", type="string", help="PDB file for the cluster_center")
parser.add_option("-n", "--clust", dest="clust", type="int", help="cluster_index")
parser.add_option("-f", "--frames", dest="frames", type="int", help="Number of frames")
parser.add_option("-s", "--start_frame", dest="start_frame", type="int", help="Number of frames")
parser.add_option("-o", "--output", dest="outfile", type="string", help="Output log file")

(options, args) = parser.parse_args()
print "Setting up GIST calculations."
lig = structure.StructureReader(options.clust_center_file).next()

#gridcntr = sum(lig.getXYZ())/len(lig.atom)
gridcntr = lig.getXYZ()[options.clust]
#print gridcntr
############################################################################
# Edit this section for different systems
gridspacn = [ 0.5, 0.5, 0.5 ]
#gridcntr = np.array([13.95, 14.59, 15.08]) # dimensions
griddim = [ 40, 40, 40 ]

############################################################################
g = Gist(options.cmsname, options.trjname, gridcntr, gridspacn, griddim)
gist_logfile = open(options.outfile + "_desmond_gist_energy.log", "w") 
gist_logfile.write("#Grid setup for the system in DX header format:\n")
gist_logfile.write('# Data calculated by the VMD volmap function\n')
gist_logfile.write('object 1 class gridpositions counts %d %d %d\n' % (g.grid.shape[0], g.grid.shape[1], g.grid.shape[2])) 
gist_logfile.write('origin %.1f %.1f %.1f\n' % (g.origin[0], g.origin[1], g.origin[2]))
gist_logfile.write('delta %.1f 0 0\n' % (g.spacing[0]))
gist_logfile.write('delta 0 %.1f 0\n' % (g.spacing[1]))
gist_logfile.write('delta 0 0 %.1f\n' % (g.spacing[2]))
gist_logfile.write('object 2 class gridconnections counts %d %d %d\n' % (g.grid.shape[0], g.grid.shape[1], g.grid.shape[2]))
gist_logfile.write('object 3 class array type double rank 0 items %d data follows\n' % (g.grid.shape[0]*g.grid.shape[1]*g.grid.shape[2]))
gist_logfile.write("#EndHeader\n")
print "Performing energy calculations ..."
t = time.time()
g.getVoxelEnergies(options.frames, options.start_frame)
print "energy calcs took seconds.", time.time() - t
g.normalizeVoxelQuantities(options.frames, gist_logfile)
g.getDipoleDotProduct()

t = time.time()
#print "Performing entropy calculations ..."
#g.getVoxelEntropies(options.frames, 0.5, gist_logfile)
#print "entropy calcs took seconds.", time.time() - t    
#g.writeGistData(options.outfile)
g.writeGistDipoleData(options.outfile)
print "Getting most probable config ..."
#g.getMostProbableConfig(options.frames, 0.5, gist_logfile)


gist_logfile.close()

# commands to run
#cd /data/Complementarity_Project_Data/systems/methane/new_gist_test_calcs
#$SCHRODINGER/run /data/Dropbox/gist-desmond-v3/gist_desmond_energy.py -i ../simulation/md_trial_04/04pd/04_pd.cms -t ../simulation/md_trial_04/04pd/j04_pd_mth_trj/ -l ../simulation/md_trial_04/04pd/methane.mae -f 10 -d 0 -o test
#gcc -O3 -lm -fPIC -shared -I /home/kamran/desmond/mmshare-v24012/lib/Linux-x86_64/include/python2.7 -I /home/kamran/desmond/mmshare-v24012/lib/Linux-x86_64/lib/python2.7/site-packages/numpy/core/include/ -o _gistcalcs.so _gistcalcs.c

# example
# $SCHRODINGER/run /data/Dropbox/desmond-gist/run_gist.py -i /microwayhome/kamran/water_structure_project_simulations/casp3/casp3_pd_2/casp3_pd_2-out.cms -t /microwayhome/kamran/water_structure_project_simulations/casp3/casp3_pd_2/casp3_pd_2_trj/ -l ../clustering/03_organizewaters/clustercenterfile.pdb -n 0 -f 10000 -s 2000 -o test
