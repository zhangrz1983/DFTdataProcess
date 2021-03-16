'''
2020.09.05
Tried:
    pymatgen PWOutput, only returns energy
    pymatgen PWInput.from_file and .from_string, both report error
Therefore, use ase

When use ase, got the following error
>>> a = read_espresso_out('vcrl.out')
>>> for i in a:
...     print(type(i))
... 
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "python3.8/site-packages/ase/io/espresso.py", line 169, in read_espresso_out
    for image_index in image_indexes:
TypeError: 'int' object is not iterable

dig the code, found that image_indexes should be int value, as index = -1
        image_indexes = all_config_indexes[index]

so change line 169 as the following, i.e. add [] to image_indexes to make it a one item
list, and it works well
    for image_index in [image_indexes]:
'''

import os
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.espresso import read_espresso_out

from CrystalToolkit.geometry import planar_structure_normalization, \
                                    change_lattice_constants, \
                                    multiple_structure_match


calc_dir = '/data/scratch/exw735/pwscf.output/'
not_planar_list = []
not_run_list = []
stru_list = {}

for dir in os.listdir(calc_dir):
    pw_out = calc_dir + dir
    strus_ase = read_espresso_out(pw_out)
    if strus_ase is not None:
        # strus_ase is a one-item list
        stru_ase = strus_ase[-1]
        #for stru_ase in strus_ase:
        stru = AseAtomsAdaptor.get_structure(stru_ase)
        stru.sort()

        is_planar, stru = planar_structure_normalization(structure=stru)
        if is_planar:
            newlatt = [[None,None,None], [None,None,None], [None,None,15.0]]
            stru = change_lattice_constants(structure=stru, lattice=newlatt)

            compo = stru.composition.reduced_composition.alphabetical_formula
            if compo not in stru_list:
                stru_list[compo] = []
            if not multiple_structure_match(structure=stru, structure_list=stru_list[compo]):
                stru.to(fmt='poscar', filename=dir+'.vasp')
                stru_list[compo].append(stru)
        else:
            not_planar_list.append(dir)
    else:
        not_run_list.append(dir)

with open ('/data/scratch/exw735/not_planar', 'w') as f:
    for item in not_planar_list:
        f.write("%s\n" % item)

with open ('/data/scratch/exw735/not_run', 'w') as f:
    for item in not_run_list:
        f.write("%s\n" % item)