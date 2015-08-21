# -*- coding: utf-8 -*-
"""
Spyder Editor
"""
from optparse import OptionParser
from schrodinger import structure
import numpy as np

_version = "$Revision: 0.0 $"
usage = """
$SCHRODINGER/run test.py input.cms input_trj ligand outfile

Description:
This script ...
"""
# parse command line options
parser = OptionParser()
parser.add_option("-l", "--gist_log_file", dest="gist_log", type="string", help="Log file output by main GIST code")
parser.add_option("-g", "--gist_data_file", dest="gist_data", type="string", help="Data file output by main GIST code")
#parser.add_option("-d", "--check_data_field", dest="check_datafield", type="int", help="column to apply cutoff")
#parser.add_option("-c", "--cutoff", dest="cutoff", type="float", help="Cutoff value for the above column")
#parser.add_option("-w", "--write_data_field", dest="dx_datafield", type="int", help="column to write in dx file")
parser.add_option("-c", "--clust_center", dest="clust_center_file", type="string", help="PDB file for the cluster_center")
parser.add_option("-n", "--clust", dest="clust", type="int", help="cluster_index")
parser.add_option("-o", "--output_prefix", dest="out_name", type="string", help="Prefix for DX file")

#cmsname = "/microwayhome/kamran/water_structure_project_simulations/abl/abl_pd_2/abl_pd_2-out.cms"
#trjname = "/microwayhome/kamran/water_structure_project_simulations/abl/abl_pd_2/abl_pd_2_trj/"
#dsim = create_simulation(cmsname, trjname)

(options, args) = parser.parse_args()

lig = structure.StructureReader(options.clust_center_file).next()
clust_cntr = lig.getXYZ()[options.clust]
clust_cntr = np.array([-9.24, -4.42, -5.32])

check_datafield = 4
dx_datafield = 15
# 7, 11, 15
cutoff = 0
# obtain head for dx files
dx_header = open(options.gist_log, "r").readlines()[1:9]
#print dx_header

# obtain column names to gnerate corresponding dx files
data = open(options.gist_data, "r").readlines()
gist_header = data[0]
data_keys = gist_header.strip("\n").split()
# obtain voxel data, where each voxel index serves as the key
voxeldata = {}
for l in data[1:]:
    float_converted_data = [float(x) for x in l.strip("\n").split()[1:]]
    voxeldata[int(l.strip("\n").split()[0])] = float_converted_data
    #print int(l.strip("\n").split()[0])
    #print float_converted_data
# specify column of the gist output whose value will be checked for cutoff
#options.check_datafield = options.check_datafield
# specify column of the gist output for which dx file will be generated
#options.dx_datafield = options.dx_datafield + 1

# go over each column of the gist output
#print data_keys
#print "Generating dx file for %s using a cutoff value of %f for %s" % (data_keys[dx_datafield+1], cutoff, data_keys[check_datafield+1]) 

f = open(options.out_name+"_"+data_keys[dx_datafield+1]+".dx",'w')
# write lines for dx header
for header_line in dx_header:
    f.write(header_line)

# access voxels three at a time
for i in range(0, len(voxeldata.keys()), 3):
    key_sublist = voxeldata.keys()[i:i+3] # create a list of voxel ids for the three voxels currently accessed 
    # for each voxel id
    for l in key_sublist:
        point = np.asarray(voxeldata[l][0:3])
        dist = np.linalg.norm(clust_cntr-point)
        #if voxeldata[l][check_datafield] >= cutoff and dist <= 3.5:
        #if voxeldata[l][check_datafield] >= cutoff and dist > 3.5 and dist <= 5.5:
        if voxeldata[l][check_datafield] >= cutoff and dist > 5.5 and dist <= 8.5:
            #print point
            f.write("%0.8f " % voxeldata[l][dx_datafield]) # access relevant value for this voxel
        else:
            f.write("%0.1f " % 0.0) # access relevant value for this voxel
    f.write("\n")
#f.write("object 'occupancy (all)' class field\n")
f.close()


# example:
# python /data/Dropbox/desmond-gist/generate_high_occ_dipole_dx.py -l desmond-gist-energy.log -g gist_data_test.txt -d 4 -c 2 -w 8 -o clust0_dipole_dp_high_occ
