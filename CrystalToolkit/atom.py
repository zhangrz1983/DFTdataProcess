
import numpy as np
from itertools import product


def bond_length_estimation(atom_list, elem_data, factor=1.0):
    '''
    Estimate the maximum distance between any two types of atoms in  
    atom_list using crystal radius of the atoms

    Args:
        atom_list: a list containing atomic symbols, like ['Na','Cl']
        elem_data: atom properties from the json file
        factor: scale factor 
    Return:
        the value of estimated maximum bond length
    '''
    pair_list = list(product(atom_list, repeat=2))
    bondlen_list = []
    for pair in pair_list:
        bondlen = 0
        for atom in pair:
            bondlen += elem_data[atom]['crystal_rad']
        bondlen_list.append(bondlen)
    return np.array(bondlen_list).max() * factor


def average_atomic_energy(atomic_formula, elem_data):
    '''
    This can be used for all the structures, not only high-entropy
    Calculated the average atomic energy of all atoms in a structure, 
    by summing all the atomic energies and then divided by number of atoms

    Args:
        atomic_formula: a dictionary, in the format of pymatgen formula
        elem_data: atom properties from the json file
    Return:
        average atomic energy in eV
    '''
    total_atoms_ene = 0.0
    n_atoms = 0
    for atom, num in atomic_formula.items():
        total_atoms_ene += (elem_data[atom]['atom_energy'] * num)
        n_atoms += num
    return total_atoms_ene / n_atoms


def atoms_rad_statis(atomic_formula, elem_data):
    '''
    Give atomic raduis statistics, i.e. standard deviation and (max-min) in a structure. 
    Currently it is designed for high-entropy ceramics materails, the statistics is on cations
    And the anion atom is found using the anion_num value

    Args:
        atomic_formula: a dictionary, in the format of pymatgen formula
        elem_data: atom properties from the json file
    Return:
        cation and anion atom list
        standard deviation and (max-min) of atomic and crystal radius 
    '''
    # assum anions have the maximum number, mostly this is true
    anion_num = max(atomic_formula.values())

    cation = []
    atom_rads = []
    cryst_rads = []
    for atom, num in atomic_formula.items():
        # the if to make sure only run atomic statistics on cations
        if num == anion_num:
            anion = atom
        else:
            # sometimes the cations is not equimolar, so use num to loop
            for i in range(int(num)):
                cation.append(atom)
                atom_rads.append(elem_data[atom]['atomic_rad'])
                cryst_rads.append(elem_data[atom]['crystal_rad'])
    return anion, cation, \
            np.array(atom_rads).std(), np.array(cryst_rads).std(), \
            np.array(atom_rads).ptp(), np.array(cryst_rads).ptp()


