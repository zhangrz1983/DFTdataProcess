
import os
import numpy as np
import math
from pymatgen import Structure, Lattice
from pymatgen.io.vasp.inputs import Poscar

from CrystalToolkit.composition import solid_solution_random
from CrystalToolkit.structure import poscar_fix_top_bottom, stacking_fault_generation


# note that slip is only along a1vect_uvw direction
slip_systems = [
    {
        'hkl': [1, 1, 0], 
        'a1vect_uvw': [1, -1, 0], 
        'a2vect_uvw': [0, 0, 1]
    },
    {
        'hkl': [1, 1, 1], 
        'a1vect_uvw': [1, -1, 0], 
        'a2vect_uvw': [1, 1, -2]
    },
    {
        'hkl': [1, 1, 1], 
        'a1vect_uvw': [1, 1, -2], 
        'a2vect_uvw': [1, -1, 0]
    },
]

anions = ['N']
cations = {
    'Ti': 4.238, 
    'Zr': 4.622, 
    'Nb': 4.444, 
    'Ta': 4.410,
#    'Ti-Zr-Nb-Ta': 4.428,
}

a2 = 0.0
for anion in anions:
    for cation in cations.keys():  
        for slip_sys in slip_systems:
            hkl, a1vect_uvw, a2vect_uvw = slip_sys.values()
            sizemults = [1, 1, 4]
            a1_max, nstep = 0.25, 11
            if a1vect_uvw == [1, 1, -2]:
                a1_max, nstep = 0.5, 21
            if len(cation) > 2:
                a1_max, nstep = 2, 21

            for a1 in np.linspace(0, a1_max, nstep):
                # every time stru should be built
                # otherwise it seems to be changed by the subroutine 
                # probably consider to use istructrue?
                if len(cation) < 3:
                    stru = Structure.from_spacegroup('Fm-3m', Lattice.cubic(cations[cation]),
                                                [cation, anion],[[0.5, 0.5, 0.5], [0, 0, 0]])
                    stru_fault = stacking_fault_generation(structure=stru, sizemults=sizemults, 
                                                    a1=a1, a2=a2, hkl=hkl, a1vect_uvw=a1vect_uvw, 
                                                    a2vect_uvw=a2vect_uvw)
                else:
                    sizemults = [x1 * x2 for x1, x2 in zip(sizemults, [2,2,1])]
                    hec_elems = cation.split('-')
                    # pymatgen seems not accept fraction in the from_spacegroup method
                    # only this can be done stru['X'] = 'X0.5Y0.5', so use hec_elems[0]
                    stru = Structure.from_spacegroup('Fm-3m', Lattice.cubic(cations[cation]),
                                                [hec_elems[0], anion],[[0.5, 0.5, 0.5], [0, 0, 0]])
                    stru_fault = stacking_fault_generation(structure=stru, sizemults=sizemults, 
                                                    a1=a1, a2=a2, hkl=hkl, a1vect_uvw=a1vect_uvw, 
                                                    a2vect_uvw=a2vect_uvw)
                    stru_fault = solid_solution_random(structure=stru_fault, 
                                                        elem_frac_site=hec_elems[0], 
                                                        elem_list=hec_elems)

                stru_fault.sort()
                selective_dynamics = poscar_fix_top_bottom(stru_fault)
                fault_sys_poscar = Poscar(structure=stru_fault, 
                                            selective_dynamics=selective_dynamics)

                compo_name = cation + '-' + anion
                if math.isclose(a1, 0.0):
                    dir = compo_name + '_' \
                            + '_'.join([''.join(map(str, hkl)), 'plane'])
                else:
                    dir = compo_name + '_' \
                            + '_'.join([''.join(map(str, hkl)), 'plane', 
                                        ''.join(map(str, a1vect_uvw))]) \
                            + '_' + format(a1, '.2f')

                if not os.path.isdir(dir):
                    os.mkdir(dir)
                    fault_sys_poscar.write_file(dir + '/POSCAR', ) 


