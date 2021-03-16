'''


Do not find primitive cell, because c lattice is considered, so in some prim_cell beta and gamma != 90
And posit is from the conv_cell, so this sometimes make error

'''
import numpy as np
import os
import shutil

from pymatgen import Structure

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from CrystalToolkit.geometry_properties import planar_structure_normalization
from CrystalToolkit.geometry_properties import change_lattice_constant, multiple_structure_match
from CrystalToolkit.mongodb_retrieval import get_mongodb_entries


elem_list = [ 'B', 'C',
              'Al','Si','P', 'S',
              'Ga','Ge','As','Se',
              'In','Sn','Sb','Te','I' ]
#elem_list = ['C']

work_dir = '/workspace/RZ/CharDen/post_tran_elem/chgden_2dnets/'
if not os.path.exists(work_dir):
    os.mkdir(work_dir)

for elem in elem_list :
    docs = get_mongodb_entries(dir_prefix='relax_2dnets/' + elem +'_')

    stru_protos = []
    for doc in docs :
        stru = Structure.from_dict(doc['output']['crystal'])
        #stru = SpacegroupAnalyzer(stru, symprec=0.1).find_primitive()
        if_planar, stru = planar_structure_normalization(structure=stru)
        if if_planar:
            stru = change_lattice_constant(structure=stru)
            duplic_stru = multiple_structure_match(structure=stru, structure_list=stru_protos)

            if not duplic_stru : 
                stru_protos.append(stru)
                calc_dir = work_dir + doc['dir_name'].split('/')[-1]
                os.mkdir(calc_dir)
                stru.to(fmt = 'poscar', filename = calc_dir + '/POSCAR')

