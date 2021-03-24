'''

'''

import os
import json
import sys
import pandas as pd
import numpy as np
from pymatgen import Structure, Lattice

from CrystalToolkit.structure import structure_from_json, elements_sorted_filename
from CrystalToolkit.composition import composition_from_one_list


struct_file = os.path.join(sys.path[0], '../CrystalToolkit/structure.json')
with open(struct_file) as f:
    structure_list = json.load(f)

stru_componts = structure_list['entropy_forming_ability']['five_cations_rocksalt_pocc']['structures']
# The lattice constants read from the file is 1/2 of the cubic crystal cell parameter
# in OQMD they are stored as ((a,a,0)(a,0,a)(0,a,a)) for the primitive cell
# therefore "a" is half the length of cubic cell parameter
lat_vals = pd.read_csv('/home/rz/Desktop/B1_NaCl_lattice_constant_OQMD.csv', index_col=0)

ntype = 5
cation_list = [
            'Al', 'Si',  
            'Ti', 'V',  'Cr',
            'Zr', 'Nb', 'Mo', 
            ]
anion_list = [
            'N', 
            ]

cation_comp_lists = composition_from_one_list(elem_list=cation_list, ntype=ntype)
for cation_comp_list in cation_comp_lists:
    # get lattice constant of the solid solution by averaging all the binary
    cell_lat = 0
    for cation in cation_comp_list:
        cell_lat += lat_vals[cation][anion_list[0]]
    cell_lat = cell_lat / len(cation_comp_list) 
    
    for i, stru_compont in enumerate(stru_componts):
        stru = structure_from_json(composition_list=list(cation_comp_list)+anion_list, 
                                structure_components=stru_compont)
        
        # re-scale the lattice constant to make the geometry optimization easier
        # note that ntype *2 ( = 10 ) is used to scale the volume
        # 'cell_lat' is half the lattic parameter of cubic NaCl-type structure
        # so the volume is only for one atom
        # while the POCC structures read from json file is for 10 atoms 
        scale_factor = ( cell_lat ** 3 * (ntype * 2) / stru.volume ) ** (1./3.)

        # in older version of pymatgen, stru.lattice cannot be change, so revise the code
        # to change structure
        #latt = stru.lattice.as_dict()
        #latt['matrix'] = (np.array(stru.lattice.as_dict()['matrix']) * scale_factor).tolist()
        #latt = Lattice.from_dict(latt)
        #stru.lattice = latt

        d_stru = stru.as_dict()
        d_stru['lattice']['matrix'] = (np.array(d_stru['lattice']['matrix']) * scale_factor).tolist()
        stru = Structure.from_dict(d_stru)

        # write into vasp poscar file using the naming for high-entropy materials
        f = elements_sorted_filename(structure=stru, incl_num=True) + '_' + str(i+1)
        os.mkdir(f)
        stru.to(fmt='poscar', filename=f+'/POSCAR')
