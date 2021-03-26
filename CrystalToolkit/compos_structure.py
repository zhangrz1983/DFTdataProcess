

import os
import random
import math
from itertools import permutations, combinations, product
from pymatgen.alchemy.materials import TransformedStructure

try:
    from pymatgen.transformations.advanced_transformations import SQSTransformation
except:
    print('Old version of pymatgen is used. SQSTransformation cannot be imported')

def solid_solution_sqs(structure, elem_frac_site, elem_frac_comp, sqs_scaling):
    '''
    Use pymatgen SQSTransformation (which call ATAT mcsqs program)
    to generate disordered structure. 

    Args:
        structure: pymatgen structure
        elem_frac_site: the factional occupy site in the original structure, e.g. 'Ti'
        elem_frac_comp: solid solution composition, e.g. 'Ti0.25Zr0.25Hf0.25Nb0.25'
        sqs_scaling (int or list): (same as pymatgen scaling in SQSTransformation)
                    Scaling factor to determine supercell. Two options are possible:
                    a. (preferred) Scales number of atoms, e.g., for a structure with 8 atoms,
                       scaling=4 would lead to a 32 atom supercell
                    b. A sequence of three scaling factors, e.g., [2, 1, 1], which
                       specifies that the supercell should have dimensions 2a x b x c
    Return:
        pymatgen structure, SQS                   
    '''
    # build another pymatgen structure
    structure[elem_frac_site] = elem_frac_comp
    ts = TransformedStructure(structure, [])

    # the directory must be set in SQSTransformation, otherwise the work dir 
    # will be changed by this function
    workdir = os.getcwd()
    ts.append_transformation(SQSTransformation(scaling=sqs_scaling, search_time=1, 
                                directory=workdir, reduction_algo=False))
    return ts.structures[-1]


def solid_solution_random(structure, elem_frac_site, elem_list):
    '''
    An alternative way to generate random structrues, where SQS has limitations.
    The method is typicall acceptable when the supercell is large, i.e. hundreds of atoms

    For example, when use for atomman stacking fault generations:
    If mcsqs is used after the fault is generated, then the atoms at the fault have different
    geometry enviorement from the bulk atoms. If mcsqs is applied to the surface system without
    fault, after the mcsqs in pymatgen the structure will past to atomman, and the surface method
    must be used again, the atomman fault method will bulid two identical SQS above and under the 
    fault (slide plane), even sizemults=[1,1,1], i.e. this double the surface system

    This function use random.shuffle(), alternatively random.choice() can be used. However, 
    random.choice (code immediately below) does not give equal number of two types of atoms, 
    sometimes the discrpency is large
        if str(fault_sys_pymatgen[i]).split()[-1] == elem_frac_site:
            fault_sys_pymatgen[i] = random.choice(elem_list)

    Args:
        structure: pymatgen structure
        elem_frac_site: the site for solid solution
        elem_list: species and composition of the solid solution
    Return:
        pymatgen supercell structure with random occupuation 
    '''
    
    # only cations have fraction occuputation, so the length of atom_list is half of the fault system
    atom_list = []
    for i in range(len(elem_list)):
        atom_list = atom_list + ( [i] * int(len(structure) / 2 / len(elem_list) + 0.5) )

    random.shuffle(atom_list)
    j = 0
    for i in range(len(structure)):
        if str(structure[i]).split()[-1] == elem_frac_site:
            structure[i] = elem_list[atom_list[j]]
            j += 1

    return structure









