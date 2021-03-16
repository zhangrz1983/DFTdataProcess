# coding: utf-8
# Distributed under the terms of the MIT License.
'''
Collection of functions generate/modify structure input files 
typically for VASP or PWscf
'''


import os
import math
from pymatgen.io.pwscf import PWInput
from pymatgen import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.espresso import read_espresso_out
try:
    import atomman
except:
    print('atomman not installed')

from CrystalToolkit.geometry import planar_structure_normalization, \
                                    change_lattice_constants


def elements_sorted_filename(structure, incl_num=False):
    '''
    Name Calculat dir or filename using the format, like
    H_He or (if incl_num=True) 01_H_02_He
    This way shows the elements and order the compositions

    Args:
        structure: pymatgen structure
        incl_num: whether include atomic number of the elements
                in the file name. Sometimes include number helps
                to sort file names. 
    Return:
        dir or file name include sorted atomic symbols, 
        e.g. 08_O_12_Mg or O_Mg
    '''
    site_dict = {}
    for site in structure.sites:
        site_dict[site.specie.number] = site.specie.name

    fname = []
    for num, symb in sorted(site_dict.items()):
        if incl_num:
            fname.append(format(num, '02d'))
        fname.append(symb)

    return '_'.join(fname)


def poscar_fix_top_bottom(structure):
    '''
    For stacking fault system, which is typicall a surface system, fixing the top and bottom
    layers during the geometry optimization.  
    Currently only z direction is allowed to relax

    Args:
        structure: pymatgen structure
        elem_frac_site: the site for solid solution
        elem_list: species and composition of the solid solution
    
    Return:
        3*N selective_dyanmics matrix as a list for pymatgen Poscar method
        (N is number of atoms in the pymatgen struture)
    '''
    
    # fix the top and bottom layer during geometry optimization
    frac_coords = structure.frac_coords
    max_z = frac_coords[:,2].max()
    min_z = frac_coords[:,2].min()
    tol = 0.001

    selective_dynamics = []
    for i in range(len(structure)):
        if abs(frac_coords[i][2] - max_z) < tol or abs(frac_coords[i][2] - min_z) < tol:
            selective_dynamics.append([False, False, False])
        else:
            selective_dynamics.append([False, False, True])
    return selective_dynamics


def pwscf_write_input(structure, calculation='vc-relax', pseudo_suffix='_pbe_gbrv1.4.UPF',
                        pseudo_dir='/data/home/exw735/Software/QE/pseudo',
                        ecutwfc=30.0, degauss=0.03, kspace=0.25):
    '''
    This function write PWscf input file using pymatgen PWInput class.  
    Required fields are add into the dictionaries

    Args:
        structure: pymatgen structure
    Return:
        pymatgen PWInput object
    '''
    control_dict = {'calculation': calculation, 'pseudo_dir': pseudo_dir}
    pseudo_dict = {}
    for atom in list(structure.composition.get_el_amt_dict().keys()):
        pseudo_dict[atom] = atom + pseudo_suffix

    system_dict = {'occupations': 'smearing', 'degauss': degauss, 
                    'ecutwfc': ecutwfc, }
    cell_dict = {'cell_factor': 8.0}
    kgrid = (
        max(1, int(2 * math.pi / (kspace * (structure.lattice.a) ) + 0.5)), 
        max(1, int(2 * math.pi / (kspace * (structure.lattice.b) ) + 0.5)), 
        max(1, int(2 * math.pi / (kspace * (structure.lattice.c) ) + 0.5))
    )
    return PWInput(structure=structure, pseudo=pseudo_dict, control=control_dict, 
                    system=system_dict, electrons=None, ions=None, cell=cell_dict, 
                    kpoints_mode='automatic', kpoints_grid=kgrid)


def pwscf_output_structure(pw_out):
    '''
    2020.09.05 tried:
        pymatgen PWOutput, only returns energy
        pymatgen PWInput.from_file and .from_string, both report error
    Therefore, use ase

    Made the following change after line 160 in  
    /data/home/exw735/.local/lib/python3.8/site-packages/ase/io/espresso.py  
    Note that as 'yield' at the end of the function, although 'return' presents,  
    the function still return a 'generator'

        if results_config_indexes:
            image_indexes = [results_config_indexes[index]]
            if indexes[_PW_CONVERG]:
                try:
                    image_indexes = [results_config_indexes[index-1]]
                except:
                    return
        else:
            return

    Args:
        pw_out: quantum espresso output file
    Return:
        pymatgen structure
    '''

    # read_espresso only returns one ase atoms object
    for stru_ase in read_espresso_out(pw_out):
        stru = AseAtomsAdaptor.get_structure(stru_ase)
        return stru.sort()


def stacking_fault_generation(structure, sizemults, a1, a2, hkl, a1vect_uvw, a2vect_uvw):
    '''
    Generate stacking fault struture using Atomman.  
    Please note the below about Atomman:  
    One layer means one 'ucell' in the defect.StackingFault, along the 'hkl' direction.  
    So 'sizemults[2]=6' means 6 ucell above the stacking fault slide plane, and 6 layers below.  
    i.e. totally 12 ucell/layers in the surface system

    Args:
        structure: pymatgen structure
        sizemults: how to extend the struture when building surface
        a1: fraction of the slide movement along a1vect_uvu
        a2: fraction of the slide movement along a2vect_uvu
        hkl: the stacking fault slide plane
        a1vect_uvw: the stacking fault moving direction
        a2vect_uvw: the stacking fault moving direction

    Return:
        pymatgen structure containing stacking fault
    '''
    
    # generate stakcing fault object using atomman, and then change to pymatgen for random occupuation
    stru_atomman = atomman.load('pymatgen_Structure', structure)
    stack_fault = atomman.defect.StackingFault(hkl=hkl, ucell=stru_atomman, 
                                                a1vect_uvw=a1vect_uvw, a2vect_uvw=a2vect_uvw)
   
    # it is a bit wired that surface_sys is not used in the rest of the code,
    # but that is how atomman is designed
    surface_sys = stack_fault.surface(shift=stack_fault.shifts[0], even=True, 
                                        sizemults=sizemults, vacuumwidth=15)
    fault_sys = stack_fault.fault(a1=a1, a2=a2)
    fault_sys_pymatgen = Structure.from_str(fault_sys.dump('poscar'), fmt='poscar')
    return fault_sys_pymatgen


def structure_from_json(composition_list, structure_components):
    '''
    Read crystal structure component in the stored json format.  
    
    Args:
        composition_list: a list elements in the structure
        structure_components: read from the json file
    Return:
        pymatgen structure
    '''
    #if len(structure_components['species']) == len(composition_list):
    species = []
    for num in structure_components['species']:
        species.append(composition_list[num-1])
    lattice = structure_components['lattice']
    coords = structure_components['coords']

    return Structure(species=species, lattice=lattice, coords=coords)

    # The following is for ASE
    # nacl = crystal(symbols=['Na','Cl'], basis=strus['NaCl'].frac_coords, cell=strus['NaCl'].lattice.matrix)

