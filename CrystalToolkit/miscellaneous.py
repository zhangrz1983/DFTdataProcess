
import os
import numpy as np
import json
from pymatgen import Structure

def poscars_to_json():
    '''
    Read the poscars in the EFA nat comm paper
    https://doi.org/10.1038/s41467-018-07160-7
    Convert into json format to put into structure.json
    '''
    structs = []
    for poscar in os.listdir('.'):
        s = Structure.from_file(poscar)
        struct = {}
        lattice = s.lattice.matrix.tolist()
        species = []
        coords = []
        for site in s.sites:
            coords.append(list(site.frac_coords))
            species.append(site.species_string)
        species = sorted(list(map(int, '-'.join(species).replace('C','6').replace('Ti','1').replace('Zr','2').replace('Hf','3').replace('Nb','4').replace('Ta','5').split('-'))))
        struct['species'] = species
        struct['lattice'] = lattice
        struct['coords'] = coords
        structs.append(struct)

    with open('../1.json', 'w') as f:
        json.dump(structs, f)


def oqmd_lattice_to_csv(infile='B1_NaCl_lattice_constant_OQMD'):
    '''
    Read the OQMD into CVS file
    '''
    periodic_table = [
        'H',  'He', 
        'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 
        'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 
        'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
        'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe', 
        'Cs', 'Ba', 
                    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
                        'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 
                    'Ac', 'Th', 'Pa', 'U',  'Np', 'Pu']

    df = pd.read_csv(infile, sep='\s+| |\t', index_col=[0,1], header=None, engine='python')
    df = df.unstack()
    df.columns = df.columns.droplevel()
    df.columns = [ x.replace('1','') for x in df.columns.values]
    df.index = [ x.replace('1','') for x in df.index.values]
    df = df.fillna(0)
    df = df.reindex(columns=periodic_table, index=periodic_table, fill_value=0) 
    df = df + df.transpose()
    df.to_csv(infile + '.csv')

