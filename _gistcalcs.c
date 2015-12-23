/*
  C helper module for intensive calcs.
#=============================================================================================
# VERSION CONTROL INFORMATION
#=============================================================================================
__version__ = "$Revision: $ $Date: $"
# $Date: $
# $Revision: $
# $LastChangedBy: $
# $HeadURL: $
# $Id: $

#=============================================================================================

#=============================================================================================
# INSTALLATION INSTRUCTIONS
#=============================================================================================

  To compile on Linux:

  gcc -O3 -lm -fPIC -shared -I(directory with Python.h) -I(directory with numpy/arrayobject.h) -o _gistcalcs.so _gistcalcs.c

  For a desmond installation of python 2.5 (change path up to desmond directory, rest should be the same):
  
  gcc -O3 -lm -fPIC -shared -I /home/kamran/desmond/mmshare-v24012/lib/Linux-x86_64/include/python2.7 -I /home/kamran/desmond/mmshare-v24012/lib/Linux-x86_64/lib/python2.7/site-packages/numpy/core/include/ -o _gistcalcs.so _gistcalcs.c
  gcc -O3 -lm -fPIC -shared -I /opt/schrodinger2014-2/mmshare-v26017/lib/Linux-x86_64/include/python2.7/ -I /opt/schrodinger2014-2/mmshare-v26017/lib/Linux-x86_64/lib/python2.7/site-packages/numpy/core/include/ -o _gistcalcs.so _gistcalcs.c
  
  For a default installation of python 2.5:
  
  gcc -O3 -lm -fPIC -shared -I/usr/local/include/python2.5 -I/usr/local/lib/python2.5/site-packages/numpy/core/include -o _gistcalcs.so _gistcalcs.c

/Users/Kamran/anaconda/include/python2.7/Python.h
/Users/Kamran/anaconda/pkgs/numpy-1.9.0-py27_0/lib/python2.7/site-packages/numpy/core/include/
*/
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double dist_mic(double* x1, double* x2, double* x3, double* y1, double* y2, double* y3, double* b1, double* b2, double* b3) {
    /* Method for obtaining inter atom distance using minimum image convention
     */
    double dx, dy, dz;
    dx = *x1 - *y1;
    dy = *x2 - *y2;
    dz = *x3 - *y3;
    if (dx > *b1/2.0) dx -= *b1; 
    else if (dx < -*b1/2.0) dx += *b1; 
    if (dy > *b2/2.0) dy -= *b2;
    else if (dy < -*b2/2.0) dy += *b2;
    if (dz > *b3/2.0) dz -= *b3; 
    else if (dz < -*b3/2.0) dz += *b3;
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

    }
    
double dist(double x1, double x2, double x3, double y1, double y2, double y3) {
    /* Method for Euclidean distance between two points
     */
    double dx, dy, dz;
    dx = x1-y1;
    dy = x2-y2;
    dz = x3-y3;
    return sqrt(pow(dx, 2)+ pow(dy, 2)+ pow(dz, 2));
    }

double calcEnergy(int at1, int at2, 
                PyObject * pos,
                PyArrayObject * chg,
                PyArrayObject * vdw,
                PyArrayObject * img){
    //printf("calculating energy...\n");
    double *b_x, *b_y, *b_z; // variables to store box info
    double *at1x, *at1y, *at1z; // variables to store first atom's coordinates
    double *at1c, *at1sig, *at1eps; // variables to store first atom's non-bonded params 
    double *at2x, *at2y, *at2z; // variables to store second atom's coordinates
    double *at2c, *at2sig, *at2eps; // variables to store second atom's non-bonded params 
    double comb_sig, comb_eps;
    double dist, dist6, dist12, aij, bij;
    double LJ, Elec, totalEnergy;

    b_x = (double *) PyArray_GETPTR1(img, 0); 
    b_y = (double *) PyArray_GETPTR1(img, 1);
    b_z = (double *) PyArray_GETPTR1(img, 2);
    //NOTE: This function use global ids as arguments for atom 1 and atom 2
    at1x = PyArray_GETPTR2(pos, at1, 0); // obtain x, y, z
    at1y = PyArray_GETPTR2(pos, at1, 1); 
    at1z = PyArray_GETPTR2(pos, at1, 2); 
    at1c = PyArray_GETPTR1(chg, at1); // charge on this atom
    at1sig = PyArray_GETPTR2(vdw, at1, 0);
    at1eps = PyArray_GETPTR2(vdw, at1, 1);
    //printf("water atom ID x y z charge: %i %5.3f %5.3f %5.3f %5.3f\n", at1, *at1x, *at1y, *at1z, *at1c);
    at2x = PyArray_GETPTR2(pos, at2, 0); // obtain x, y, z
    at2y = PyArray_GETPTR2(pos, at2, 1); 
    at2z = PyArray_GETPTR2(pos, at2, 2); 
    at2c = PyArray_GETPTR1(chg, at2); // charge on this atom
    at2sig = PyArray_GETPTR2(vdw, at2, 0);
    at2eps = PyArray_GETPTR2(vdw, at2, 1);
    //printf("solute atom ID x y z charge sig eps: %i %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", at2, *at2x, *at2y, *at2z, *at2c, *at2sig, *at2eps);
    // we first calculate water molecule's oxygen non-bonded interactions (both vdw and elec) with every solute atom
    dist = dist_mic(at1x, at1y, at1z, at2x, at2y, at2z, b_x, b_y, b_z); // distance (based on minimum image convention)
    // Coulombic interaction calculation
    Elec = ((*at1c)*(*at2c))/dist; 
    //printf("Energy: %f\n", *(double *)PyArray_GETPTR2(voxel_data, v_id, 12));
    // Coulombic interaction calculation
    comb_sig = (*at1sig + *at2sig)/2;
    comb_eps = sqrt((*at1eps)*(*at2eps));
    dist6 = pow(dist, 6);
    dist12 = dist6 * dist6;
    aij = 4*comb_eps*pow(comb_sig, 12);
    bij = -4*comb_eps*pow(comb_sig, 6);
    LJ = (aij/dist12)+(bij/dist6);
    totalEnergy = Elec + LJ;
    //printf("Energy between %i and %i: %f\n", at1, at2, totalEnergy);
    //Py_DECREF(pos);
    /*
    Py_DECREF(b_x);
    Py_DECREF(b_y);
    Py_DECREF(b_z);
    Py_DECREF(at1x);
    Py_DECREF(at1y);
    Py_DECREF(at1z);
    Py_DECREF(at1c);
    Py_DECREF(at1sig);
    Py_DECREF(at1eps);
    Py_DECREF(at2x);
    Py_DECREF(at2y);
    Py_DECREF(at2z);
    Py_DECREF(at2c);
    Py_DECREF(at2sig);
    Py_DECREF(at2eps);
    */
    return totalEnergy;
    }

void energy( double * ox, double * oy, double * oz,
                PyArrayObject * target_atom_ids,
                PyArrayObject * wat_O_ids,
                PyObject * all_coords,
                PyArrayObject * all_charges,
                PyArrayObject * all_vdw,
                PyArrayObject * img_info,
                PyArrayObject * voxel_data,
                PyArrayObject * wat_info,
                int * w_id, int v_id){

    // energy calculation is done in a separate function but for oxygen-oxygen pair, distance and energy are calculated within this
    // function because of requirements for Enbr calculations
    double *wx, *wy, *wz; // variables to store water atom's coordinates (only for H-atoms)
    double *w_c, *w_sig, *w_eps; // variables to store target water oxygen atom's non-bonded params (separately for O and H-atoms)
    double *o_c, *o_sig, *o_eps; // variables to store current water oxygen atom's non-bonded params (separately for O and H-atoms)
    double comb_sig, comb_eps;
    double dist, dist6, dist12, aij, bij;

    double *b_x, *b_y, *b_z; // variables to store box info

    // following variables store info that helps retrieving pseudo interaction sites
    int n_atomic_sites, n_pseudo_sites, wat_begin_id, pseudo_begin_id, oxygen_index, offset, offset2;
    // following variables control looping over atoms
    int n_solute_atoms, n_wat_O_atoms;
    int wat_at_i, wat_ps_i, solute_at_i, solvent_at_i, solvent_ps_i;
    int wat_h, solvent_h;
    int * target_at;
    // initialize energy values for this water molecule
    double E_sw = 0.0;
    double E_ww = 0.0;
    double E_nbr = 0.0;
    int nbrs = 0;
    b_x = (double *) PyArray_GETPTR1(img_info, 0); 
    b_y = (double *) PyArray_GETPTR1(img_info, 1);
    b_z = (double *) PyArray_GETPTR1(img_info, 2);

    n_atomic_sites = *(int *)PyArray_GETPTR1(wat_info, 0);
    n_pseudo_sites = *(int *)PyArray_GETPTR1(wat_info, 1);
    wat_begin_id = *(int *)PyArray_GETPTR1(wat_info, 2);
    pseudo_begin_id = *(int *)PyArray_GETPTR1(wat_info, 3);
    oxygen_index = *(int *)PyArray_GETPTR1(wat_info, 4);
    offset = (*w_id - 1 - wat_begin_id)/n_atomic_sites;
    //printf("Number of atom sites: %i\n", n_atomic_sites);
    //printf("Number of pseudo sites: %i\n", n_pseudo_sites);
    //printf("Water atom site start ID: %i\n", wat_begin_id);
    //printf("Water pseudo start site ID: %i\n", pseudo_begin_id);
    //printf("General oxygen index: %i\n", oxygen_index);

    n_solute_atoms = PyArray_DIM(target_atom_ids, 0);
    n_wat_O_atoms = PyArray_DIM(wat_O_ids, 0);
    // Begin E_sw calculations
    // iterate over all solute atoms
    for (solute_at_i = 0; solute_at_i < n_solute_atoms; solute_at_i++){
        // iterate over each atomic interaction site for this water
        //printf("Solute atom with atomic index %i\n", solute_at_i);
        for (wat_at_i = 0; wat_at_i < n_atomic_sites; wat_at_i++){
            //printf("Water atom site %i with atomic index %i and global index %i\n", wat_at_i, w_id, w_id - 1); // gid = w_id - 1 + wat_at_i - oxygen_index
            E_sw += calcEnergy(*w_id - 1 + wat_at_i, solute_at_i, all_coords, all_charges, all_vdw, img_info);
            }
        // iterate over each pseudo interaction site for this water
        for (wat_ps_i = 0; wat_ps_i < n_pseudo_sites; wat_ps_i++){
            //printf("Water pseudo site %i with global index %i\n", wat_ps_i, pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i); // gid = w_id - 1 + wat_at_i - oxygen_index
            offset = (*w_id - 1 - wat_begin_id)/n_atomic_sites;
            E_sw += calcEnergy(pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i, solute_at_i, all_coords, all_charges, all_vdw, img_info);
            }
        }
    *(double *)PyArray_GETPTR2(voxel_data, v_id, 12) += E_sw;
    //printf("E_sw: %f\n", E_sw);

    // End E_sw calculations

    // Begin E_wat calculations
    // begin iterating over each solvent oxygen atom
    for (solvent_at_i = 0; solvent_at_i < n_wat_O_atoms ; solvent_at_i++){

        double E_wat = 0.0;

        // 1. Calculate distance and interaction energy of current water O-atom and target water O-atom
        target_at = (int *) PyArray_GETPTR1(wat_O_ids, solvent_at_i); // get atom id
        if (*target_at-1 == *w_id-1) continue; // if global ids match, skip calculation
        wx = PyArray_GETPTR2(all_coords, *target_at-1, 0); // obtain x, y, z using global id
        wy = PyArray_GETPTR2(all_coords, *target_at-1, 1); 
        wz = PyArray_GETPTR2(all_coords, *target_at-1, 2); 
        w_c = PyArray_GETPTR1(all_charges, *target_at-1); // charge on this atom
        w_sig = PyArray_GETPTR2(all_vdw, *target_at-1, 0);
        w_eps = PyArray_GETPTR2(all_vdw, *target_at-1, 1);
        // obtain non-bonded params for the current water O atom
        o_c = PyArray_GETPTR1(all_charges, *w_id-1); // charge on this atom
        o_sig = PyArray_GETPTR2(all_vdw, *w_id-1, 0);
        o_eps = PyArray_GETPTR2(all_vdw, *w_id-1, 1);

        //printf("target water O atom ID x y z charge sig eps: %i %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", *target_at-1, *wx, *wy, *wz, *w_c, *w_sig, *w_eps);
        //printf("current water O atom ID x y z charge sig eps: %i %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n", w_id-1, *ox, *oy, *oz, *o_c, *o_sig, *o_eps);
        
        dist = dist_mic(ox, oy, oz, wx, wy, wz, b_x, b_y, b_z); // distance (based on minimum image convention)
        // Coulombic interaction calculation
        E_wat += ((*o_c)*(*w_c))/dist; 
        //printf("electrostatics: %f for charges %f and %f at distance %f\n", ((*o_c)*(*w_c))/dist, *o_c, *w_c);
        comb_sig = (*o_sig + *w_sig)/2;
        comb_eps = sqrt((*o_eps)*(*w_eps));
        dist6 = pow(dist, 6);
        dist12 = dist6 * dist6;
        aij = 4*comb_eps*pow(comb_sig, 12);
        bij = -4*comb_eps*pow(comb_sig, 6);
        E_wat += (aij/dist12)+(bij/dist6);
        //printf("LJ: %f\n", ((aij/dist12)+(bij/dist6)));
        //printf("Energy between %i and %i: %f\n", w_id-1, *target_at-1, ((*o_c)*(*w_c))/dist + (aij/dist12)+(bij/dist6));

        // 2. Iterate over H-atoms of current water to calculate energy between target water's O-atom & current water's H-atoms
        // begin iterating over current water's hydrogen atoms
        for (wat_h = *w_id; wat_h <= *w_id+1; wat_h++){
            //printf("target water's H atom ID: %i\n", solvent_h);
            E_wat += calcEnergy(wat_h, *target_at-1, all_coords, all_charges, all_vdw, img_info);
            } // end iterating over current water's hydrogen atoms
        // 3. Iterate over target water's H-atoms to calculate interactions with current water's O-atom, H-atoms and pseudo-sites
        // begin iterating over target water's hydrogen atoms
        for (solvent_h = *target_at; solvent_h <= *target_at+1; solvent_h++){
            //printf("target water's H atom ID: %i\n", solvent_h);
            E_wat += calcEnergy(solvent_h, *w_id-1, all_coords, all_charges, all_vdw, img_info);
            // begin iterating over current water's hydrogen atoms
            for (wat_h = *w_id; wat_h <= *w_id+1; wat_h++){
                //printf("current water's H atom ID: %i\n", wat_h);
                E_wat += calcEnergy(solvent_h, wat_h, all_coords, all_charges, all_vdw, img_info);
                } // end iterating over current water's hydrogen atoms
            // begin iterating over each pseudo interaction site for current water
            for (wat_ps_i = 0; wat_ps_i < n_pseudo_sites; wat_ps_i++){
                //printf("Water pseudo site %i with global index %i\n", wat_ps_i, pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i); // gid = w_id - 1 + wat_at_i - oxygen_index
                E_wat += calcEnergy(solvent_h, pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i, all_coords, all_charges, all_vdw, img_info);
                } // end iterating over each pseudo interaction site for current water

            } // end iterating over target water's hydrogen atoms
        // 4. Iterate over current water's pseudo sites to calculate interaction with target water's O-atom
        // with this step current water's pseudo sites' interactions with all atomic sites of target water are finished.
        // begin iterating over each pseudo interaction site for current water 
        for (wat_ps_i = 0; wat_ps_i < n_pseudo_sites; wat_ps_i++){
            //printf("Current Water pseudo site %i with global index %i\n", wat_ps_i, pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i); // gid = w_id - 1 + wat_at_i - oxygen_index
            E_wat += calcEnergy(pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i, *target_at-1, all_coords, all_charges, all_vdw, img_info);
            } // end iterating over each pseudo interaction site for current water
        
        // 5. Iterate over target water's pseudo sites to calculate interactions with current water's O-atom, H-atoms and pseudo sites
        // begin iterating over each pseudo interaction site for target water
        offset2 = (*target_at - 1 - wat_begin_id)/n_atomic_sites;
    
        for (solvent_ps_i = 0; solvent_ps_i < n_pseudo_sites; solvent_ps_i++){
            //printf("Target Water pseudo site %i with global index %i\n", solvent_ps_i, pseudo_begin_id + offset*n_pseudo_sites + solvent_ps_i); // gid = w_id - 1 + wat_at_i - oxygen_index
            E_wat += calcEnergy(pseudo_begin_id + offset2*n_pseudo_sites + solvent_ps_i, *w_id-1, all_coords, all_charges, all_vdw, img_info);
            // begin iterating over current water's hydrogen atoms
            for (wat_h = *w_id; wat_h <= *w_id+1; wat_h++){
                //printf("current water's H atom ID: %i\n", wat_h);
                E_wat += calcEnergy(wat_h, pseudo_begin_id + offset2*n_pseudo_sites + solvent_ps_i, all_coords, all_charges, all_vdw, img_info);
                } // end iterating over current water's hydrogen atoms
            // *** insert pseudo site loop here *** 
            // begin iterating over each pseudo interaction site for current water
            for (wat_ps_i = 0; wat_ps_i < n_pseudo_sites; wat_ps_i++){
                //printf("Water pseudo site %i with global index %i\n", wat_ps_i, pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i); // gid = w_id - 1 + wat_at_i - oxygen_index
                E_wat += calcEnergy(pseudo_begin_id + offset2*n_pseudo_sites + solvent_ps_i, pseudo_begin_id + offset*n_pseudo_sites + wat_ps_i, all_coords, all_charges, all_vdw, img_info);
                } // end iterating over each pseudo interaction site for current water
            } // end iterating over each pseudo interaction site for target water
        

        // E_wat is energy of current water molecule with the target water molecule
        //printf("E_wat: %f\n", E_wat);
        E_ww += E_wat;

        if (dist <= 3.5){
            nbrs += 1;
            E_nbr += E_wat;        
            }

        } // end iterating over each solvent oxygen atom

        
    //printf("Water %i, E_sw: %f, E_ww: %f, E_nbr_tot: %f, E_nbr_norm: %f, nbrs: %i\n", w_id-1, E_sw, E_ww/2.0, E_nbr/2.0, E_nbr/(nbrs*2.0), nbrs);

    *(double *)PyArray_GETPTR2(voxel_data, v_id, 14) += E_ww;
    // End E_ww calculations

    // sanity check on Enbr, if no nbrs found, E_nbr should be zero    
    if (isnan(E_nbr/(nbrs*2.0))){
        //printf("Water %i: E_tot: %f E_shell_tot: %f E_shell_norm: %f Nbrs: %i\n", w_id, E_tot/2.0, E_shell_tot/2.0, E_shell_tot/(wat_nbr*2.0), wat_nbr);
        *(double *)PyArray_GETPTR2(voxel_data, v_id, 16) += 0;
        *(double *)PyArray_GETPTR2(voxel_data, v_id, 18) += nbrs;
        }
    else {
        *(double *)PyArray_GETPTR2(voxel_data, v_id, 16) += E_nbr/(nbrs*2.0);
        *(double *)PyArray_GETPTR2(voxel_data, v_id, 18) += nbrs;
        }
    
    Py_DECREF(ox);
    Py_DECREF(o_c);
    Py_DECREF(o_sig);
    Py_DECREF(o_eps);
    /*
    Py_DECREF(wy);
    Py_DECREF(wz);
    Py_DECREF(w_c);
    Py_DECREF(w_sig);
    Py_DECREF(w_eps);
    Py_DECREF(b_x);
    Py_DECREF(b_y);
    Py_DECREF(b_z);
    */
    } // end energy calculation




PyObject *_gistcalcs_processGrid(PyObject *self, PyObject *args)
{
    // variables reterived from python
    int frames, start_frame, num_atoms; // total number of frames and atoms in the system
    PyObject *sendWatCoords, *arglist_sendWatCoords, *call_result;// python function, its argument tuple and variable to store its results
    // Following variables are pointers to arrays that come from Python and contain various important pieces of info
    PyObject *coords; // coordinates for all atoms, retrieved through a python callback, are stored in this array
    PyArrayObject *all_at_ids, *solute_at_ids, *wat_oxygen_ids, *wat_all_ids;
    PyArrayObject *charges, *vdw, *box, *grid_dim, *grid_orig;
    PyArrayObject *voxel_data; // array for voxel, arrives here empty
    PyArrayObject *wat_index_info;
    PyArrayObject *test_array;

    // variables to be used locally
    int i_frames, i_wat; // frame counter, water counter
    int n_wat;    
    double grid_max_x, grid_max_y, grid_max_z;
    double grid_orig_x, grid_orig_y, grid_orig_z;
    int grid_dim_x, grid_dim_y, grid_dim_z;
    double grid_index_x, grid_index_y, grid_index_z; 
    double *wat_x, *wat_y, *wat_z;
    double *h1x, *h1y, *h1z, *h2x, *h2y, *h2z;
    int *wat_id;
    double wat_dist;
    int n_atomic_sites, n_pseudo_sites, wat_begin_id, pseudo_begin_id, oxygen_index;
    int voxel_id;

    double wat_translated_x, wat_translated_y, wat_translated_z;

    // Argument parsing to reterive everything sent from Python correctly    
    if (!PyArg_ParseTuple(args, "iOO!O!O!O!O!O!O!O!O!O!O!O!:processGrid",
                            &num_atoms, &sendWatCoords,
                            &PyArray_Type, &coords,
                            &PyArray_Type, &wat_index_info,
                            &PyArray_Type, &all_at_ids,
                            &PyArray_Type, &solute_at_ids,
                            &PyArray_Type, &wat_oxygen_ids,
                            &PyArray_Type, &wat_all_ids,
                            &PyArray_Type, &charges,
                            &PyArray_Type, &vdw,
                            &PyArray_Type, &box,
                            &PyArray_Type, &grid_dim,
                            &PyArray_Type, &grid_orig,
                            &PyArray_Type, &voxel_data))
        {
            return NULL; /* raise argument parsing exception*/
        }
    // Consistency checks for Python Callback

    if (!PyCallable_Check(sendWatCoords))
    {
        PyErr_Format(PyExc_TypeError,
                     "function is not callcable");
        return NULL;
    }

    // Set grid max points and grid origin
    grid_max_x = *(double *)PyArray_GETPTR1(grid_dim, 0) * 0.5 + 1.5;
    grid_max_y = *(double *)PyArray_GETPTR1(grid_dim, 1) * 0.5 + 1.5;
    grid_max_z = *(double *)PyArray_GETPTR1(grid_dim, 2) * 0.5 + 1.5;
    grid_orig_x = *(double *)PyArray_GETPTR1(grid_orig, 0);
    grid_orig_y = *(double *)PyArray_GETPTR1(grid_orig, 1);
    grid_orig_z = *(double *)PyArray_GETPTR1(grid_orig, 2);
    grid_dim_x = (int)*(double *)PyArray_GETPTR1(grid_dim, 0);
    grid_dim_y = (int)*(double *)PyArray_GETPTR1(grid_dim, 1);
    grid_dim_z = (int)*(double *)PyArray_GETPTR1(grid_dim, 2);
    //printf("grid origin: %f %f %f \n", grid_orig_x , grid_orig_y, grid_orig_z);
    //printf("grid max: %f %f %f \n", grid_max_x , grid_max_y, grid_max_z);
    //printf("grid dim: %i %i %i \n", (int)*(double *)PyArray_GETPTR1(grid_dim, 0), (int)*(double *)PyArray_GETPTR1(grid_dim, 1), (int)*(double *)PyArray_GETPTR1(grid_dim, 2));

    // Parse index information array
    n_atomic_sites = *(int *)PyArray_GETPTR1(wat_index_info, 0);
    n_pseudo_sites = *(int *)PyArray_GETPTR1(wat_index_info, 1);
    wat_begin_id = *(int *)PyArray_GETPTR1(wat_index_info, 2);
    pseudo_begin_id = *(int *)PyArray_GETPTR1(wat_index_info, 3);
    oxygen_index = *(int *)PyArray_GETPTR1(wat_index_info, 4);
    n_wat = PyArray_DIM(wat_oxygen_ids, 0);
    //printf("nwat: %i \n", n_wat);

    // voxel_id initialized to zero (for each water it's voxel ID will be stored in this variable) 
    // Now we need to iterate over every water atom inside grid for its voxel assignment
    for (i_wat = 0; i_wat < n_wat; i_wat ++) 
    {
        wat_id = (int *) PyArray_GETPTR1(wat_oxygen_ids, i_wat); // obtain index for this atom (this is not array index, this is unique atom id)
        // use water ID to get the correct x, y, z coordinates from coord array
        wat_x = (double *)PyArray_GETPTR2(coords, *wat_id-1, 0);
        wat_y = (double *)PyArray_GETPTR2(coords, *wat_id-1, 1); 
        wat_z = (double *)PyArray_GETPTR2(coords, *wat_id-1, 2);
        wat_translated_x = *wat_x - grid_orig_x;
        wat_translated_y = *wat_y - grid_orig_y;
        wat_translated_z = *wat_z - grid_orig_z;
        //printf("water oxygen ID %i and coordinates %f %f %f\n", *wat_id, *wat_x, *wat_y, *wat_z);
        // check if the distance between wateer coordinates and grid origin is less than the max grid point
        // this means do calculations only waters inside the grid
        if (wat_translated_x <= grid_max_x && wat_translated_y <= grid_max_y && wat_translated_z <= grid_max_z &&
            wat_translated_x >= -1.5 && wat_translated_y >= -1.5 && wat_translated_z >= -1.5)
        {
            if (wat_translated_x >= 0 && wat_translated_y >= 0 && wat_translated_z >= 0)
            {
                // transform water coordinates in units of grid dimensions
                grid_index_x = (int) (wat_translated_x/0.5);
                grid_index_y = (int) (wat_translated_y/0.5);
                grid_index_z = (int) (wat_translated_z/0.5);
                // check if water coords (in grid dimensions) are less than grid dimensions in each direction
                if (grid_index_x < grid_dim_x && grid_index_y < grid_dim_y && grid_index_z < grid_dim_z)
                {
                    // obtain the voxel ID for this water
                    //voxel_id = ((int)grid_index_x*(int)*(double *)PyArray_GETPTR1(grid_dim, 1) + (int)grid_index_y)*(int)*(double *)PyArray_GETPTR1(grid_dim, 2) + (int)grid_index_z;
                    voxel_id = (grid_index_x*grid_dim_y + grid_index_y)*grid_dim_z + grid_index_z;
                    // Energy calculations
                    //energy_ww(wat_x, wat_y, wat_z, wat_oxygen_ids, coords, charges, vdw, box, voxel_data, *wat_id, voxel_id, n_wat);
                    // get hydrogen atom coords
                    h1x = (double *)PyArray_GETPTR2(coords, *wat_id, 0);
                    h1y = (double *)PyArray_GETPTR2(coords, *wat_id, 1); 
                    h1z = (double *)PyArray_GETPTR2(coords, *wat_id, 2);
                    h2x = (double *)PyArray_GETPTR2(coords, *wat_id + 1, 0);
                    h2y = (double *)PyArray_GETPTR2(coords, *wat_id + 1, 1); 
                    h2z = (double *)PyArray_GETPTR2(coords, *wat_id + 1, 2);
                    //printf("Energy value for this voxel: %f\n",*(double *)PyArray_GETPTR2(voxel_data, voxel_id, 13));
                    //printf("water coords %f %f %f\n", *wat_x, *wat_y, *wat_z);
                    // send water x, y, z to python
                    arglist_sendWatCoords = Py_BuildValue("(iddddddddd)", voxel_id, *wat_x, *wat_y, *wat_z, *h1x, *h1y, *h1z, *h2x, *h2y, *h2z);
                    call_result = PyEval_CallObject(sendWatCoords, arglist_sendWatCoords);
                    Py_DECREF(arglist_sendWatCoords);
                    Py_DECREF(call_result);
                    energy(wat_x, wat_y, wat_z, solute_at_ids, wat_oxygen_ids, coords, charges, vdw, box, voxel_data, wat_index_info, wat_id, voxel_id);
                   
                    //printf("water coords %f %f %f\n", *wat_x, *wat_y, *wat_z);
                    //printf("grid indices %f %f %f\n", grid_index_x, grid_index_y, grid_index_z);
                    //printf("grid indices %i %i %i\n", (int)grid_index_x, (int)grid_index_y, (int)grid_index_z);
                    //printf("(%i*%i + %i)*%i + %i = ", (int)grid_index_x, (int)*(double *)PyArray_GETPTR1(grid_dim, 1), (int)grid_index_y, (int)*(double *)PyArray_GETPTR1(grid_dim, 2), (int)grid_index_z);
                    //printf("voxel id: %i\n", voxel_id);
                    // Once voxel id is obtained, it is used to reterieve various properties of the voxel and modify them based
                    // based on the calculations. Here, we reterive voxel water population and raise it by 1
                    *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 4) += 1.0;
                    // insert code here that assigns hydrogens correctly
                    *(double *)PyArray_GETPTR2(voxel_data, voxel_id, 6) += 2.0;
                    //*(int *) PyArray_GETPTR2(voxel_data, voxel_id, 3) += 1.0;
                    // 
                    //printf("voxel id: %i\n", voxel_id);
                }
            }
        }
        //printf("wat_id x y z max_x max_y max_z: %i %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n", *wat_id, *wat_x, *wat_y, *wat_z, *grid_max_x, *grid_max_y, *grid_max_z);
    }
    Py_DECREF(coords);
    Py_DECREF(sendWatCoords);
    
    /*
    Py_DECREF(wat_x);
    Py_DECREF(wat_y);
    Py_DECREF(wat_z);
    Py_DECREF(h1x);
    Py_DECREF(h1y);
    Py_DECREF(h1z);
    Py_DECREF(h2x);
    Py_DECREF(h2y);
    Py_DECREF(h2z);
    */
    //Py_DECREF(all_at_ids);
    //Py_DECREF(solute_at_ids);
    //Py_DECREF(wat_oxygen_ids);
    //Py_DECREF(wat_all_ids);
    //Py_DECREF(charges);
    //Py_DECREF(vdw);
    //Py_DECREF(box);
    //Py_DECREF(grid_dim);
    //Py_DECREF(grid_orig);
    //Py_DECREF(voxel_data);
    //Py_DECREF(wat_index_info);

    return Py_BuildValue("i", 1);
}

PyObject *_gistcalcs_getNNOr(PyObject *self, PyObject *args){
    int nwtot, n, l;
    double NNor, totR, dW, wat_or_ent;
    double voxel_dTSor = 0.0;
    double rx, ry, rz;
    PyArrayObject *voxel_wat_Eulers; 
    double twopi = 2*M_PI;

    // Argument parsing to reterive everything sent from Python correctly    
    if (!PyArg_ParseTuple(args, "iO!:getNNOr",
                            &nwtot,
                            &PyArray_Type, &voxel_wat_Eulers))
        {
            return NULL; /* raise argument parsing exception*/
        }
    // for each water in the voxel
    for (n = 0; n < nwtot; n++){
        NNor = 10000;
        for (l = 0; l < nwtot; l++){
            if(l == n) continue;
            //printf("Calculating orientational distancce between water: %i and %i\n", l, n);
            rx = cos(*(double *)PyArray_GETPTR2(voxel_wat_Eulers, l, 0)) - cos(*(double *)PyArray_GETPTR2(voxel_wat_Eulers, n, 0));
            ry = *(double *)PyArray_GETPTR2(voxel_wat_Eulers, l, 1) - *(double *)PyArray_GETPTR2(voxel_wat_Eulers, n, 1);
            rz = *(double *)PyArray_GETPTR2(voxel_wat_Eulers, l, 2) - *(double *)PyArray_GETPTR2(voxel_wat_Eulers, n, 2);
            if      (ry>M_PI) ry = twopi-ry;
            else if (ry<-M_PI) ry = twopi+ry;
            if      (rz>M_PI) rz = twopi-rz;
            else if (rz<-M_PI) rz = twopi+rz;
            dW = sqrt(rx*rx + ry*ry + rz*rz);
            //dR = 0.0;
            // get six-D distance
            // get translational nearest neighbor            
            if (dW>0 && dW<NNor) NNor = dW;
            // get six-D nearest neighbor            
        }
        //calculate translational entropy
        if (NNor<9999 && NNor>0) {
            //printf("Nearest neighbour translational distance: %f\n", NNtr);
            //wat_tr_ent = log(nwtot*NNtr*NNtr*NNtr/(3.0*twopi));
            //voxel_dTStr_norm += wat_tr_ent;
            wat_or_ent = log(nwtot*NNor*NNor*NNor/(3.0*twopi));
            voxel_dTSor += wat_or_ent;
            
        }

        
    }
    // 
    //*(double *)PyArray_GETPTR1(ent, 2) += voxel_dTSor_norm;
    return Py_BuildValue("f", voxel_dTSor);

}


/* Method Table
 * Registering all the functions that will be called from Python
 */

static PyMethodDef _gistcalcs_methods[] = {
    {
        "processGrid",
        (PyCFunction)_gistcalcs_processGrid,
        METH_VARARGS,
        "Process grid"
    },
    {
        "getNNOr",
        (PyCFunction)_gistcalcs_getNNOr,
        METH_VARARGS,
        "get voxel orientational entropy"
    },
    {NULL, NULL}
};

/* Initialization function for this module
 */

PyMODINIT_FUNC init_gistcalcs() // init function has the same name as module, except with init prefix
{
    // we produce name of the module, method table and a doc string
    Py_InitModule3("_gistcalcs", _gistcalcs_methods, "Process GIST calcs.\n");
    import_array(); // required for Numpy initialization
}
