
import os
import sys
import json
import argparse
import numpy as np
from pymatgen import Structure

from CrystalToolkit.atom import average_atomic_energy, atoms_rad_statis
from CrystalToolkit.geometry import bond_length_statis, check_coords
from pyMongoDB.mongodb import get_mongodb_entries


def group_data_use_composition(docs):
    '''
    When the docs have been read from one dir, this dir may have different compositions
    with different configurations (atomic configuations of supercells). 
    This function put all the configurations with one composition together. 
    
    Args:
        docs:
    Return:
        all_composits: a dictionary, key is the pretty formula, and the values
                are the mongodb documents for that composition/pretty formula
        not_finish_calc: 
    '''
    all_composits = {}
    not_finish_calc = []
    for doc in docs:
        if doc['state'] == 'killed' or doc['state'] == 'unsuccessful':
            not_finish_calc.append(doc['dir_name'].split('/')[-1])
        else:
            # use pretty formula as the key of the list
            # as different compositions have different pretty formula
            if doc['pretty_formula'] in all_composits:
                all_composits[doc['pretty_formula']].append(doc)
            else:
                all_composits[doc['pretty_formula']] = []
                all_composits[doc['pretty_formula']].append(doc)
    return all_composits, not_finish_calc        


def check_coords_for_compostions(all_composits, coord_num):
    '''

    Args:

    Return:

    '''
    # remove the structures no longer be the prototype
    # i.e. somehow amphorous
    new_all_composits = {}
    not_keep_struct = []
    for pretty_formula, docs in all_composits.items():
        for doc in docs:
            struct = Structure.from_dict(doc['output']['crystal'])
            if check_coords(structure=struct, coord_num=coord_num):
                if pretty_formula in new_all_composits:
                    new_all_composits[pretty_formula].append(doc)
                else:
                    new_all_composits[pretty_formula] = []
                    new_all_composits[pretty_formula].append(doc)
            else:
                not_keep_struct.append(doc['dir_name'].split('/')[-1])
    return new_all_composits, not_keep_struct


def get_values_for_each_composition(all_composits, properties=['final_energy_per_atom']):
    '''

    
    Args:
        all_composits: mongodb documents grouped using composition, 
                the output of function group_data_use_composition
        properties: a list of properties of interest, can be the keys in mongodb (matgendb output)
                or the properties implenment below:
                    bond_length_statis: statistics of all the bond length values
    Return:
        value_list: values of all the supercell configurations for a given composition
                the format is like
                    {
                        'TiZrNbTaC4':{                  # the composition / pretty formula
                            'final_energy_per_atom':[   # the property
                                1.234                   # value of each supercell configuration
                                2.345
                                ......
                            ]
                            'bond_length_statis':[      # the property
                                (1.2, 2.3, 3.4)         # the tuple of mean, std and ptp
                                (4.5, 5.6, 6.7)
                                ......
                            ]
                        }
                    }
    '''
    all_values = {}
    for pretty_formula, docs in all_composits.items():
        all_values[pretty_formula] = {}
        # creat empty list for each property
        for struct_property in properties:
            all_values[pretty_formula][struct_property] = []
            for doc in docs:
                if struct_property in doc['output']:
                    va = doc['output'][struct_property]
                elif struct_property in doc['analysis']:
                    va = doc['analysis'][struct_property]
                # some properties not included in the mongodb goes here
                elif struct_property == 'bond_length_statis':
                    struct = Structure.from_dict(doc['output']['crystal'])
                    va = bond_length_statis(structure=struct)
                else:
                    print('unknown property')

                all_values[pretty_formula][struct_property].append(va)
            
    return all_values


def entropy_forming_ability(all_values):
    '''
    Calculate entropy forming ability using the stardard deviation of energies 
    of a series of supercell structures of a give composition. See:
        Sarker, P., et al. (2018). Nat Commun 9(1): 4980.
        "High-entropy high-hardness metal carbides discovered by entropy descriptors." 

    The average formation enthalpy is also calculated to give an estimation of the
    Gibbs free energy in a plot. See:
        Wang, Y., et al. (2020). Advanced Theory and Simulations 3(9): 2000111.
        "Enhanced Hardness in High-Entropy Carbides through Atomic Randomness." 

    Args:
        all_values: a dictionary wit keys of pretty formula and values of the interesting
                values of all the supercell configurations of the corrersponding
                pretty formula (i.e. composition)
                output of the function get_values_of_one_composition
    Return:
        efa_values: a diction have the following format
            {
                'pretty formula':[  # 3 values for every compostion
                    [
                        number of supercell used in the calcuation
                        entropy forming ability (ene_list.std)
                        average total energy (ene_list.mean, not formation enthalpy)                        
                    ], 
                    ......
                ]
            }
    '''
    efa_values = {}
    for pretty_formula, struct_property in all_values.items():
        # get the list of energies of all supercell configurations
        all_energies = np.array(struct_property['final_energy_per_atom'])
        # get efa and average total energy
        efa_values[pretty_formula] = [len(all_energies), 
                                        1.0 / all_energies.std(), 
                                        all_energies.mean()]
    return efa_values


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='entorpy forming ability',
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--dir_prefix', type=str, default='nitride_hec5_relax',
                        #'HighEntro/entropy_forming_ability/sulfide_zinc_blende/'
                        help='part of the dir as the criterion to search mongodb')
    parser.add_argument('--coord_num', type=int, default=4,
                        help='coordinate number in the structure, default is 4 for zinc blende')

    args = parser.parse_args()
    dir_prefix = args.dir_prefix
    # this value is used to check if the geometry optimzation breaks the prototype structure
    # by find if all the atoms keep the original coordinates
    # see function check_coords for details
    coord_num = args.coord_num

    # read element related data from file
    # this file has some information that pymatgen does not have
    with open(os.path.join(sys.path[0], '../CrystalToolkit/atom.json')) as f:
        cryst_atom_data = json.load(f)

    docs = get_mongodb_entries(dir_prefix=dir_prefix)

    # group the documents with the same composition
    all_composits, not_finish_calc = group_data_use_composition(docs)
    with open('not_finish_calc', 'w') as f:
        print(not_finish_calc, file=f)

    # remove the relaxed structure becoming somehow amorphous
    all_composits, not_keep_struct = check_coords_for_compostions(all_composits, coord_num=coord_num)
    with open('not_keep_struct', 'w') as f:
        print(not_keep_struct, file=f)

    # get energy per atom and other propertiers composition-wise
    all_values = get_values_for_each_composition(all_composits, 
            properties=['final_energy_per_atom', 'bond_length_statis'])
            # properties=['max_force', 'delta_volume','percent_delta_volume'] 
    efa_values = entropy_forming_ability(all_values)

    for pretty_formula in efa_values.keys():
        # the docs for each composition have the same 'reduce_cell_formula'
        # so only take the first one
        doc = all_composits[pretty_formula][0]
        # sum over all the atomic energies
        avg_atoms_ene = average_atomic_energy(atomic_formula=doc['reduced_cell_formula'], 
                                    elem_data=cryst_atom_data['elem'])
        # calcuation formation energy
        average_energy_pa = efa_values[pretty_formula][2]
        formation_ene_pa = average_energy_pa - avg_atoms_ene         
        
        # get element related statistics, so only formula (not structure) is needed
        anion, cation, \
        cation_atom_rad_std, cation_cryst_rad_std, \
        cation_atom_rad_ptp, cation_cryst_rad_ptp \
            = atoms_rad_statis(atomic_formula=doc['reduced_cell_formula'], 
                                elem_data=cryst_atom_data['elem'])    
        cation = ','.join(sorted(cation))

        # get bond length analysis
        bond_len_vals = np.array(all_values[pretty_formula]['bond_length_statis'])
        # the mean at the end is the average over all the supercell configurations
        bond_len_mean = bond_len_vals[:,0].mean()
        bond_len_std = bond_len_vals[:,1].mean()
        bond_len_ptp = bond_len_vals[:,2].mean()
            
        print(anion, cation, efa_values[pretty_formula][0],
            round(cation_atom_rad_std, 3), round(cation_cryst_rad_std, 3),
            round(cation_atom_rad_ptp, 3), round(cation_cryst_rad_ptp, 3),
            round(formation_ene_pa, 3), round(efa_values[pretty_formula][1], 3),
            round((bond_len_mean), 3), round(bond_len_std, 3), round(bond_len_ptp, 3),
            sep=', ', flush=True)

 