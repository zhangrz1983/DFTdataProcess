'''
usage:
python 2dnets_pt_row3.py `echo $HOME`
'''

import os
import time
import json
import sys
from itertools import product
from CrystalToolkit.geometry import planar_structure_normalization, change_lattice_constants, \
                                    multiple_structure_match, minimum_bond_length
from CrystalToolkit.structure import pwscf_write_input, structure_from_json, elements_sorted_filename
from CrystalToolkit.atom import bond_length_estimation
from CrystalToolkit.composition import composition_from_two_lists, fit_noequi_sites, \
                                        composition_from_one_list


home_dir = sys.argv[1]
with open(home_dir + '/DFTdataProcess/CrystalToolkit/atom.json') as f:
    cryst_atom_data = json.load(f)

with open(home_dir + '/DFTdataProcess/CrystalToolkit/structure.json') as f:
    structure_list = json.load(f)

cation_list = [
            'Li', 'Be', 
            'Na', 'Mg', 'Al', 
            'K',  'Ca', 'Ga', 
            'Rb', 'Sr', 'In', 'Sn', 
            'Cs', 'Ba' 
            ]
anion_list = [
            'B', 'C',  'N',  'O',  'F', 
                 'Si', 'P',  'S',  'Cl',
                 'Ge', 'As', 'Se', 'Br'
            ]
pt_row_3 = ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']

for net_name, stru_compont in structure_list['2dnets'].items():
    noequi_sites = max(stru_compont['species'])
    if noequi_sites > 1 and noequi_sites < 4:
        #elem_lists = composition_from_two_lists(elem_list_1=cation_list, elem_list_2=anion_list)
        elem_lists = composition_from_one_list(elem_list=pt_row_3, ntype=2)
        for elem_list in elem_lists:
            stru_list = []
            comp_list = fit_noequi_sites(elem_list=elem_list, noequi_sites=noequi_sites)
            for j, atoms in enumerate(comp_list):
                stru = structure_from_json(composition_list=atoms, structure_components=stru_compont)
                stru.sort()
                if not multiple_structure_match(structure=stru, structure_list=stru_list):
                    # rescale the structure to aviod too close atoms
                    min_bondlen = minimum_bond_length(structure=stru)
                    atom_list = list(stru.composition.get_el_amt_dict().keys())
                    est_bondlen = bond_length_estimation(atom_list=atom_list, 
                                                elem=cryst_atom_data['elem'], factor=0.7)
                    if min_bondlen < est_bondlen:
                        scale_factor = est_bondlen / min_bondlen
                        newlatt =  [[scale_factor,scale_factor,None], 
                                    [scale_factor,scale_factor,None], 
                                    [None,None,None]]
                        stru = change_lattice_constants(structure=stru, lattice=newlatt, 
                                                        scale=True)
                    
                    stru_list.append(stru)
                    # write files
                    
                    f = net_name + '_' \
                        + elements_sorted_filename(structure=stru) \
                        + '_' + str(j+1)
                    os.mkdir(f)
                    #stru.to(fmt='poscar', filename=f+'/POSCAR')
                    #stru.to(fmt='poscar', filename=f+'.vasp')

                    pw = pwscf_write_input(structure=stru)
                    pw.write_file(filename=f+'/vcrl.in')