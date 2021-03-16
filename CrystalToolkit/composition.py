
import os
import random
import math
from itertools import permutations, combinations, product
from pymatgen.alchemy.materials import TransformedStructure

try:
    from pymatgen.transformations.advanced_transformations import SQSTransformation
except:
    print('Old version of pymatgen is used. SQSTransformation cannot be imported')


def composition_from_one_list(elem_list, ntype, perc=1.0):
    '''
    General all/partial combinations of a given element list

    Args:
        elem_list: a list of elements, e.g. ['Sc','Ti',...]
        ntype: number of types of elements in the formula, e.g. for high-entropy 
                materials typically ntype = 5
        perc: percentage of the total combinations to return
    Return:
        combination of the elements
    '''
    comb_list = list(combinations(elem_list, ntype))
    if not math.isclose(perc, 1.0):
        random.shuffle(comb_list)
        num = int(len(comb_list) * perc)
        comb_list = comb_list[:num]
    return comb_list


def composition_from_two_lists(elem_list_1, elem_list_2, ntype=2, perc=1.0):
    '''
    General all/partial combinations of two given element lists.  
    Typically one list is cation and the other is anion.

    Args:
        elem_list_1: a list of cation, e.g. ['Sc','Ti',...]
        elem_list_2: a list of anion, e.g. ['O','F',...]
        ntype: number of types of elements in the formula, e.g. for high-entropy 
                materials typically ntype = 5
        perc: percentage of the total combinations to return
    Return:
        combination of the elements
    '''
    comp_list = []
    for cation in elem_list_1:
        for anion in elem_list_2:
            comp_list.append([cation, anion])
    return comp_list


def fit_noequi_sites(elem_list, noequi_sites):
    '''
    Fit n type of element into m ( m >= n ) inequivalent site in a structure.  
    All elements in elem_list should be included, so use 'set' function

    Args:
        elem_list: element list, length is n, e.g. [cation,anion]
        noequi_sites: number of inequivalent sites, m
    Return:
    '''
    comp_list = list(product(elem_list, repeat=noequi_sites))
    # [:] is to make a copy otherwise remove will not work properly
    for comp in comp_list[:]:
        if len(set(comp)) < len(elem_list):
            comp_list.remove(comp)
    return comp_list


def split_equi_sites(elem_list, noequi_sites):
    '''
    Fit n type of element into m ( m < n ) inequivalent site in a structure.  
    Currently only 'm = n - 1' is supported.  
    For two elements occupy one site, typical rock-salt like pattern is considered.

    Args:
        elem_list: element list, length is n, e.g. [cation,anion]
        noequi_sites: number of inequivalent sites, m
    Return:
    '''
    # make a 2x2x1 supercell of the structure

    # fill 0,2 sites with one element and 1,3 sites with the other


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