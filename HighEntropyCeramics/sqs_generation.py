
import numpy as np
import os
from pymatgen import Structure
from CrystalToolkit.composition_properties import solid_solution_sqs

if __name__ == '__main__':
    stru_file = 'TiN_mp-492_conventional_standard.cif'
    elem_frac_site = 'Ti' 
    sqs_scaling = [2, 2, 2]

    # sqs must be disorder structure, so delete Ti=1.0 and Al=1.0 by endpoint=False 
    # and np.delete first item, respectively
    for frac in np.delete(np.linspace(0, 1, 8, endpoint=False), 0):
        # every time stru should be read, otherwise it seems to be changed by the subroutine
        stru = Structure.from_file(stru_file)
        elem_frac_comp = 'Ti' + str(frac) + 'Al' + str(1 - frac)
        stru_sqs = solid_solution_sqs(structure=stru, elem_frac_site=elem_frac_site, 
                                        elem_frac_comp=elem_frac_comp, sqs_scaling=sqs_scaling)
        dir = 'sqs_' + elem_frac_comp
        os.mkdir(dir)        
        stru_sqs.to(fmt='poscar', filename=dir+'/POSCAR')