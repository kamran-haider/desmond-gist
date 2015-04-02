from __future__ import division
__doc__='''

#===============================================================================
#
#          FILE:  Main classes and functions implementing hydration site identification
#                 using a grid-based algorithm 
#         USAGE:  a tester script will be provided as an example 
# 
#   DESCRIPTION:  
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  Desmond, Schrodinger Python API
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:  Kamran Haider
#   Contibutors:  
#     COPYRIGHT:  
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  02/24/2014 01:07:57 PM EDT
#      REVISION:  ---
#===============================================================================

'''
_version = "$Revision: 1.0 $"
# import shrodinger modules
from schrodinger import structure
from schrodinger.trajectory.desmondsimulation import create_simulation
from schrodinger.trajectory.atomselection import select_component
from schrodinger.trajectory.atomselection import FrameAslSelection as FAS
from schrodinger.application.desmond.cms import Vdw
import schrodinger.application.desmond.ffiostructure as ffiostructure
#import schrodinger.infra.mmlist as mmlist
#import schrodinger.infra.mm as mm
# import C-helper module
import _gistcalcs as gistcalcs

# import other python modules
import numpy as np
#from scipy.spatial import KDTree, cKDTree
import time
#import math

#################################################################################################################
# Main GIST class                                                                                               #
#################################################################################################################



class Gist:
#*********************************************************************************************#
    # Initializer function
    def __init__(self, input_cmsname, input_trjname, center, resolution, dimensions):
        """
        Data members
        """
        self.cmsname = input_cmsname
        self.dsim = create_simulation(input_cmsname, input_trjname)
        self._indexGenerator()
        self._initializeGrid(center, resolution, dimensions)
        self.voxeldata, self.voxeldict = self._initializeVoxelDict()
        self.chg, self.vdw = self._getNonbondedParams()
        self.box = self._initializePBC()

#*********************************************************************************************#
    # index generator function
    def _indexGenerator(self):

        frame = self.dsim.getFrame(0)
        # atom and global indices of all atoms in the system
        self.all_atom_ids = np.arange(len(self.dsim.cst.atom))+1
        self.all_atom_gids = np.arange(len(self.dsim.cst.atom)+self.dsim.cst._pseudo_total)
        # obtain oxygen atom and global indices for all water molecules
        oxygen_sel = FAS('atom.ele O')
        all_oxygen_atoms = oxygen_sel.getAtomIndices(frame)
        water_sel = select_component(self.dsim.cst, ['solvent'])
        solvent_atoms = water_sel.getAtomIndices(frame)
        solvent_oxygen_atoms = list(set(solvent_atoms).intersection(set(all_oxygen_atoms)))
        solvent_oxygen_atoms.sort()
        self.wat_oxygen_atom_ids = np.array(solvent_oxygen_atoms, dtype=np.int)
        self.wat_oxygen_atom_gids = self.wat_oxygen_atom_ids - 1
        # obtain atom and global indices for all water atoms
        wat_id_list = self.getWaterIndices(self.wat_oxygen_atom_ids)
        self.wat_atom_ids = wat_id_list[0]
        self.wat_atom_gids = wat_id_list[1]        
        # obtain all non-water atom and global indices
        self.non_water_atom_ids = np.setxor1d(self.all_atom_ids, self.wat_atom_ids).astype(int)
        self.non_water_gids = np.setxor1d(self.all_atom_gids,self.wat_atom_gids)
 
#*********************************************************************************************#
    # retrieve water atom indices for selected oxygen indices 
    def getWaterIndices(self, oxygen_atids):
        # here we will get data that is required to create previous index mapper object
        # first step is to obtain solvent forcefield structure
        solvent_ffst = None
        # obtain solvent fsst by iterating over all 'types' of forcefield (i.e., solute, solvent, ion)
        for ffst in self.dsim.cst.ffsts:
            if ffst.parent_structure.property['s_ffio_ct_type'] == 'solvent':
                if solvent_ffst is not None:
                    raise Exception("does not support multiple solvent ct.")
                solvent_ffst = ffst

        # set oxygen index to none
        oxygen_index = None
        # set types of pseudo particles to 0
        npseudo_sites = 0
        # set number of solvent atom types to 0
        natom_sites = 0
       # for each forcefield site (which is any 'site' on the structure to which params are assigned)
        for i, site in enumerate(solvent_ffst.ffsite):
            # check if this site belongs to Oxygen atoms
            if site.vdwtype.upper().startswith('O'):
                # if oxygen index is already defined, raise exception otherwise set oxygen index to this site
                if oxygen_index is not None:
                    raise Exception("water molecule has more than two oxygen atoms")
                oxygen_index = i
            # check if this site belongs to pseudoparticle, if yes raise corresponding number
            if site.type.lower() == 'pseudo':
                npseudo_sites += 1
            # check if this site belongs to an atom, if yes raise corresponding number
            elif site.type.lower() == 'atom':
                natom_sites += 1
        # at the end of this loop we have checked all possible forcefield sites to get the correst index for oxygen
        # in addition we get total number of atoms and pseudopartciles on a solvent site (water in this case) 
        if oxygen_index is None:
            raise Exception("can not locate oxygen atom.")
        if natom_sites == 0:
            raise Exception("number of atoms is zero.")
        # here we totall number of atoms in solvent 
        nmols = len(solvent_ffst.parent_structure.atom)/natom_sites
        #print oxygen_index
        # this is atid for the first oxygen atom in water oxygen atom array
        wat_begin_atid = oxygen_atids[0]
        # gid in this case is atid - 1
        wat_begin_gid = wat_begin_atid - 1
        oxygen_gids = oxygen_atids - 1
        pseudo_begin_gid = wat_begin_gid + natom_sites*nmols
        id_list = []
        #return atids of atoms of selected water molecules.
        water_atids = []
        for oxygen_atid in oxygen_atids:
            for i in range(natom_sites):
                atid = oxygen_atid + i - oxygen_index
                water_atids.append(atid)
        id_list.append(np.array(water_atids))
        #return gids of particles (including pseudo sites) of selected water molecules.
        water_gids = []
        for oxygen_gid in oxygen_gids:
            for i in range(natom_sites):
                gid = oxygen_gid + i - oxygen_index
                water_gids.append(gid)
            # pseudo atoms are placed right after real atoms
            offset = (oxygen_gid - wat_begin_gid) / natom_sites
            for i in range(npseudo_sites):
                gid = pseudo_begin_gid + offset*npseudo_sites + i
                water_gids.append(gid)
        water_gids.sort()
        id_list.append(np.array(water_gids, dtype=np.int))
        self.oxygen_index = oxygen_index
        self.n_atom_sites = natom_sites
        self.n_pseudo_sites = npseudo_sites
        self.wat_begin_gid = wat_begin_gid
        self.pseudo_begin_gid = pseudo_begin_gid
        return id_list
#*********************************************************************************************#
    def _getNonbondedParams(self):
        # obtain LJ and Elec params
        #*********************************************************************************#
        vdw = [] # combined list of all vdw params from all ct's
        chg = [] # combined list of all charges from all ct's
        ct_list = [e for e in ffiostructure.CMSReader(self.cmsname)]
        struct_ct = ct_list[1:] # this means this works only on cms files with separate CT blocks
        # get total number of solute atoms
        for ct in struct_ct:
            #print ct.atom_total, len(ct.ffio.site), len(ct.ffio.pseudo)
            #else: # do the normal parsing
            ct_vdw   = [] # list of vdw objects for each ct
            ct_chg = []
            ct_pseudo_vdw = []
            ct_pseudo_chg = []
            n_atomic_sites = 0
            n_pseudo_sites = 0
            vdw_type = {} # dict of vdw types, Vdw object in this list are uninitialized
            for e in ct.ffio.vdwtype :
                vdw_type[e.name] = Vdw( (e.name,), e.funct, (e.c1, e.c2,) )
                #print (e.name,), e.funct, (e.c1, e.c2,)

            for e in ct.ffio.site: # for each site (i.e., an atom in most cases)
                #print e.type.lower()
                if e.type.lower() == 'pseudo':
                    ct_pseudo_vdw.append( vdw_type[e.vdwtype] ) # add to vdw list for this ct
                    ct_pseudo_chg.append(e.charge)
                    n_pseudo_sites += 1
                else:
                    ct_vdw.append( vdw_type[e.vdwtype] ) # add to vdw list for this ct
                    ct_chg.append(e.charge)
                    n_atomic_sites += 1

                #print e.index, e.charge
		        # check if this site belongs to pseudoparticle, if yes raise corresponding number
            #print n_atomic_sites, n_pseudo_sites
            ct_vdw *= int(ct.atom_total/n_atomic_sites)
            ct_chg *= int(ct.atom_total/n_atomic_sites)
            vdw.extend( ct_vdw )
            chg.extend( ct_chg)
            if n_pseudo_sites != 0:
                ct_pseudo_vdw *= int(len(ct.ffio.pseudo)/n_pseudo_sites)
                ct_pseudo_chg *= int(len(ct.ffio.pseudo)/n_pseudo_sites)
                #print int(ct.atom_total / len( ct.ffio.site ))
                vdw.extend( ct_pseudo_vdw )
                chg.extend( ct_pseudo_chg)

        chg = np.asarray(chg)*18.2223
        vdw_params = []
        for v in vdw:
            vdw_params.extend([v.c])
        vdw_params = np.asarray(vdw_params)
        
        #tmp_wat_ids = self.getWaterIndices(np.asarray([self.wat_atom_ids[0]]))[0]
        #tmp_wat_gids = self.getWaterIndices(np.asarray([self.wat_atom_ids[0]]))[1]
        #print tmp_wat_ids, tmp_wat_gids
        #print chg[tmp_wat_gids]
        #print vdw_params[tmp_wat_gids]

        #print vdw_params[self.water_atom_ids-1]
        return (chg, vdw_params)
        
#*********************************************************************************************#
    def _initializeGrid(self, center=[0.0,0.0,0.0], resolution=[1.0,1.0,1.0], dimensions=[1, 1, 1]):
        # set grid center, res and dimension
        self.center = np.array(center,dtype=np.float_)
        self.dims = np.array(dimensions)
        self.spacing = np.array(resolution,dtype=np.float_)
        o = self.center - 0.5*self.dims*self.spacing
        self.origin = np.around(o, decimals=2)
        # set grid size (in terms of total points alog each axis)
        length = np.array(self.dims/self.spacing, dtype=np.float_)
        self.grid_size = np.ceil(length / self.spacing + 1.0)
        self.grid_size = np.cast['uint32'](self.grid_size)
        # Finally allocate the space for the grid
        self.grid = np.zeros(self.dims,dtype=np.float_)


#*********************************************************************************************#
    def _initializeVoxelDict(self):
        voxel_dict = {}
        v_count = 0
        voxel_array = np.zeros((self.grid.size, 21), dtype="float64")
        #print voxel_dict_new.shape
        for index, value in np.ndenumerate(self.grid): 
            #point = grid.pointForIndex(index) # get cartesian coords for the grid point
            _index = np.array(index, dtype=np.int32)
            #point = self.spacing * _index + self._origin
            point = _index*self.spacing + self.origin + 0.5*self.spacing
            voxel_array[v_count, 1] = point[0]
            voxel_array[v_count, 2] = point[1]
            voxel_array[v_count, 3] = point[2]
            voxel_array[v_count, 0] = v_count
            
            #print voxel_dict_new[v_count, 0], voxel_dict_new[v_count, 1], voxel_dict_new[v_count, 2]
            #voxel_dict[v_count] = [[]] # create a dictionary key-value pair with voxel index as key and it's coords as
            #voxel_dict[v_count].append(np.zeros(14, dtype="float64"))
            v_count += 1
        return voxel_array, voxel_dict

#*********************************************************************************************#
    def _initializePBC(self):
        # for minimum image convention
        box_vectors = self.dsim.getFrame(0).box
        if box_vectors[0] == 0.0 or box_vectors[4] == 0.0 or box_vectors[8] == 0.0:
            print "Warning: Periodic Boundary Conditions unspecified!"
        else:
            box = np.asarray([box_vectors[0], box_vectors[4], box_vectors[8]])
        return box
#*********************************************************************************************#
    def sendCoords(self, i):
        #print "Processing frame: ", i+1, "..."
        coords = []
        # get frame structure, position array
        frame = self.dsim.getFrame(i)
        #frame_st = self.dsim.getFrameStructure(i)
        pos = frame.position
        pos_pseudo = frame.pseudo_position
        coords.extend(pos)
        if pos_pseudo.shape[0] != 0:
            coords.extend(pos_pseudo)
        coords = np.asarray(coords)
        return coords
#*********************************************************************************************#
    def getvoxelWatCoords(self, voxel_id, ox, oy, oz, h1x, h1y, h1z, h2x, h2y, h2z):
        pi = np.pi
        twopi = 2*np.pi
        wat_coords = [ox, oy, oz, h1x, h1y, h1z, h2x, h2y, h2z]
        # define the lab frame of reference
        xlab = np.asarray([1.0, 0.0, 0.0], dtype="float64")
        #ylab = np.asarray([0.0, 1.0, 0.0], dtype="float64")
        zlab = np.asarray([0.0, 0.0, 1.0], dtype="float64")
        # create array for water oxygen atom coords, and append to this voxel's coord list
        owat = np.asarray(wat_coords[0:3], dtype="float64")
        # create array for water's hydrogen 1 and 2
        h1wat = np.asarray(wat_coords[3:6], dtype="float64") - owat
        h2wat = np.asarray(wat_coords[6:9], dtype="float64") - owat
        #print owat, h1wat, h2wat
        # define water molecule's frame
        # H1 is water's x-axis, should be normalized
        xwat = np.copy(h1wat)
        xwat *= 1/(np.linalg.norm(h1wat))
        # z-axis is the cross-product of H1 and H2
        zwat = np.cross(xwat, h2wat)
        zwat *= 1/(np.linalg.norm(zwat))
        # y-axis is just perpendicular to z- and x-axis
        ywat = np.cross(zwat, xwat)
        ywat *= 1/(np.linalg.norm(ywat))
        # Now we proceed to Euler angle calculations between water and lab frame
        # we start with theta and we will use cosine formula for the dot product`
        dp = np.dot(zlab, zwat)
        # first we get theta which is angle between z-axes of two frames
        theta = np.arccos(dp)
        phi = 0
        psi = 0
        # if theta is between 0 and pi
        if theta > 1E-5 and theta < pi-1E-5:
            # define a new vector which is perpendicular to both z-axes
            node = np.cross(zlab, zwat)
            node *= 1/(np.linalg.norm(node))
            # get angle phi which is the angle between node and xlab
            dp = np.dot(node, xlab)
            if dp <= -1.0:
                phi = pi
            elif dp >= 1.0:
                phi = pi
            else:
                phi = np.arccos(dp)
                # check if angle phi is between 0 and 2pi
            if phi > 0.0 and phi < twopi:
                # define new vector v which is perpendicular to xlab and node
                v = np.cross(xlab, node)
                dp = np.dot(v, zlab)
                if dp < 0:
                    phi = twopi - phi
                    # get angle psi
                dp = np.dot(xwat, node)
            if dp <= - 1.0:
                psi = pi
            elif dp > 1.0:
                psi = 0.0
            else:
                psi = np.arccos(dp)
            if psi > 0.0 and psi < twopi:
                v = np.cross(node, xwat)
                dp = np.dot(v, zwat)
                if dp < 0:
                    psi = twopi - psi
            if not theta <= pi and theta >= 0 and phi <= twopi and phi >= 0 and psi <= twopi and psi >= 0:
                print "Error: Euler angles don't fall into range!"
                # add angle information for this water 
                #angwat = np.asarray([theta, phi, psi], dtype="float64")
        #angle_array.append([theta, phi, psi])
        # Finish translational and rotational coords for waters in this voxel
        # perform nearest neighbor search for each water for this voxel
        # first in translational space

        self.voxeldict[voxel_id][0].append([theta, phi, psi])
        #angle_array = np.asarray(angle_array, dtype="float64")
    
#*********************************************************************************************#
    # Energy Function
    def getVoxelEnergies(self, n_frame, s_frame):
        # testing new GIST module
        wat_index_info = np.array([self.n_atom_sites, self.n_pseudo_sites, self.wat_begin_gid, self.pseudo_begin_gid, self.oxygen_index], dtype="int")
        gistcalcs.processGrid(n_frame, s_frame, len(self.all_atom_ids), self.sendCoords, self.getvoxelWatCoords, 
                    wat_index_info, self.all_atom_ids, self.non_water_atom_ids, self.wat_oxygen_atom_ids, self.wat_atom_ids, self.chg, self.vdw, self.box,
                    self.dims.astype("float"), self.origin, self.voxeldata)

#*********************************************************************************************#
    def getVoxelEntropies(self, n_frame, res, logfile):
        # set constants
        gas_kcal = 0.0019872041
        rho_bulk = 0.0334
        voxel_vol = res**3
        # initialize variables to store entropy results
        dTStr_tot = 0.0
        dTSor_tot = 0.0
        total_mem = 0
        print "Memory taken up for by rest of the grid in Kbs",  sys.getsizeof(self.voxeldata)/1000
        # begin iterating over voxels
        for k in self.voxeldict.keys():
            # only for voxels with more than one water
            if self.voxeldata[k, 4] > 1.0:
                dTStr_dens = -gas_kcal*300*rho_bulk*self.voxeldata[k, 5]*np.log(self.voxeldata[k, 5])
                # voxel translational entropy density
                self.voxeldata[k, 7] = dTStr_dens
                # voxel translational entropy normalized
                self.voxeldata[k, 8] = dTStr_dens*n_frame*voxel_vol/(1.0*self.voxeldata[k, 4])
                # initialize array for storing translational and rotational coords for waters in this voxel
                angle_array = np.asarray(self.voxeldict[k][0], dtype="float64")
                total_mem += angle_array.nbytes
                dTSor_dens = gistcalcs.getNNOr(int(self.voxeldata[k, 4]), angle_array)
                # voxel orientational entropy normalized
                self.voxeldata[k, 10] = gas_kcal*300*((dTSor_dens/self.voxeldata[k, 4]) + 0.5772)
                # voxel orientational entropy density
                self.voxeldata[k, 9] = self.voxeldata[k, 10]*(self.voxeldata[k, 4]/(n_frame*voxel_vol))

                dTStr_tot += self.voxeldata[k, 7]
                dTSor_tot += self.voxeldata[k, 9]

        dTStr_tot *= voxel_vol
        dTSor_tot *= voxel_vol
        logfile.write("Total Solute-Water Translational Entropy over the grid: %.8f\n" % dTStr_tot)
        logfile.write("Total Solute-Water Orientational Entropy over the grid: %.8f\n" % dTSor_tot)
        print "Total Solute-Water Orientational Entropy over the grid: %.8f" % dTSor_tot
        print "Total Solute-Water Translational Entropy over the grid: %.8f" % dTStr_tot
        print "Memory taken up by water euler angles Kbs", total_mem/1000
#*********************************************************************************************#

    def normalizeVoxelQuantities(self, n_frame, logfile):
        voxel_vol = 0.5**3
        rho_bulk = 0.0329
        bulkwaterpervoxel = voxel_vol*rho_bulk*n_frame
        Eswtot = 0
        Ewwtot = 0
        Ewwnbrtot = 0
        nwat = 0
        for k in self.voxeldata:
            if k[4] >= 1.0:
                
                # number density of oxygen centers
                k[5] = k[4]/bulkwaterpervoxel
                # number density of hydrogens
                k[6] = k[6]/(bulkwaterpervoxel*2)
                # energy density-weighted and normalized
                k[11] = k[12]/(n_frame*voxel_vol*2.0) # E_sw_dens
                Eswtot += k[11]
                k[12] = k[12]/k[4] # E_sw_norm
                k[13] = k[14]/(n_frame*voxel_vol*2.0) # E_ww_dens
                Ewwtot += k[13]
                k[14] = k[14]/(k[4]*2.0) # E_ww_norm
                # neighbours
                k[17] = k[18]/(n_frame*voxel_vol) # nbr_dens
                k[18] = k[18]/k[4] # nbr_norm
                # E_nbr
                k[15] = k[16]/(n_frame*voxel_vol)
                k[16] = k[16]/k[4]
                Ewwnbrtot += k[15]
                enclosure = 1 - (k[18]/5.25)
                if enclosure < 0.0:
                    enclosure = 0.0
                k[19] += enclosure
                k[20] = k[19]/(n_frame*voxel_vol) # Enclosure-dens
                nwat += k[4]/(n_frame*voxel_vol)

                
        Eswtot *= voxel_vol
        Ewwtot *= voxel_vol
        Ewwnbrtot *= voxel_vol
        nwat *= voxel_vol
        print "Total number of water molecules over the grid: ", nwat
        print "Total Solute-Water Energy over the grid: ", Eswtot
        print "Total Water-Water Energy over the grid: ", Ewwtot
        print "Total Water-Water Nbr Energy over the grid: ", Ewwnbrtot
        #logfile.write("Total Solute-Water Energy over the grid: %.8f\n" % Eswtot)
        #logfile.write("Total Water-Water Energy over the grid: %.8f\n" % Ewwtot)
#*********************************************************************************************#
                    
    def writeGistData(self, outfile):
        f = open("energy_data_"+outfile+".txt", "w")
        header_2 = "voxel x y z wat g_O g_H Esw-dens Esw-norm Eww-dens Eww-norm Eww-nbr-dens Eww-nbr-norm nbr-dens nbr-norm enc-dens enc-norm\n"
        f.write(header_2)
        for k in self.voxeldata:
            l = "%i %.2f %.2f %.2f %i %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n" % \
                (k[0], k[1], k[2], k[3], k[4], k[5], k[6],
                k[11], k[12], k[13], k[14], k[15], k[16], k[17], k[18], k[19], k[20])
                #print l
            f.write(l)
        f.close()
        
