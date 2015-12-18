import numpy as np
from optparse import OptionParser
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
parser.add_option("-o", "--output_prefix", dest="out_name", type="string", help="Prefix for DX file")
(options, args) = parser.parse_args()
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

# obtain header for dx files
dx_header = open(options.gist_log, "r").readlines()[1:8]

dx_file_objects = []
#print dx_header
"""
for data_field, title in enumerate(data_keys):
    if data_field > 3:# and data_field < 19:
        print data_field, title
      
        f = open(options.out_name+"_"+title+".dx",'w')
        for header_line in dx_header:
            f.write(header_line)
        for i in range(0, len(voxeldata.keys()), 3):
            key_sublist = voxeldata.keys()[i:i+3] 
            for l in key_sublist:
                f.write("%0.6f " % (voxeldata[l][data_field-1]))
                #print v_dict[l][1][2],
            #print
            f.write("\n")
        #f.write("object 'occupancy (all)' class field\n")
        f.close()
"""
for data_field, title in enumerate(data_keys):
    #if data_field > 4:# and data_field < 6:
    #print "Writing dx file for: ", title
    if data_field > 4:
        f = open(options.out_name+"_"+title+".dx",'w')
        for line in dx_header:
            f.write(line)
        dx_file_objects.append(f)
    else:
        dx_file_objects.append(None)

#for column_i in xrange(5, len(data_keys)):
#    print column_i, dx_file_objects[column_i]
#    print voxeldata[0][column_i]


for k in xrange(0, len(voxeldata)):
    #print "writing data for voxel: ", k
    if voxeldata[k][4] > 1.0:
        for column_i in xrange(5, len(data_keys) - 1):
            #print column_i, dx_file_objects[column_i], voxeldata[k][column_i]
            dx_file_objects[column_i].write("%0.6f " % (voxeldata[k][column_i]))
            if k % 3 == 0:
                dx_file_objects[column_i].write("\n")
    else:
        for column_i in xrange(5, len(data_keys) - 1):
            dx_file_objects[column_i].write("%i " % (voxeldata[k][column_i]))
            if k % 3 == 0:
                dx_file_objects[column_i].write("\n")

for f in dx_file_objects:
    if f != None:
        f.close()
