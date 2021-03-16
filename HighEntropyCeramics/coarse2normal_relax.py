
import os
import json
import pymatgen 

from CrystalToolkit.atom import total_atomic_energy, atoms_rad_statis
from CrystalToolkit.geometry import bond_length_statis, check_coords
from CrystalToolkit.mongodb import get_mongodb_entries


def ideal_coord_num_from_dir(dir_name):
    if 'B1_NaCl' in dir_name :
        ncoord = 6
        stru_type = 'NaCl'
    if 'B2_CsCl' in dir_name :
        ncoord = 8
        stru_type = 'CsCl'
    if 'B3_ZnS' in dir_name :
        ncoord = 4
        stru_type = 'ZnS'
    return ncoord, stru_type


if __name__ == '__main__':

    with open('data.json') as json_file:
        cryst_atom_data = json.load(json_file)

    anion_num = 5.0
    docs = get_mongodb_entries(dir_prefix='HighEntro/bond_length_relax/B1_NaCl/hec5_20sqs_coarse')
    calc_dir = '/workspace/RZ/HighEntro/bond_length_relax/B1_NaCl/hec5_20sqs/'
    for doc in docs:
        total_atoms_ene = total_atomic_energy(atomic_formula=doc['reduced_cell_formula'].items(), 
                                    elem=cryst_atom_data['elem'])
        anion, cation, \
        cation_atom_rad_std, cation_cryst_rad_std, \
        cation_atom_rad_ptp, cation_cryst_rad_ptp \
            = atoms_rad_statis(atomic_formula=doc['reduced_cell_formula'].items(), 
                                    anion_num=anion_num, elem=cryst_atom_data['elem'])    
        cation = ','.join(cation)
        coord_num, stru_type = ideal_coord_num_from_dir(dir_name=doc['dir_name']) 
        if doc['state'] == 'killed' or doc['state'] == 'unsuccessful':
            print(stru_type, anion, cation, 
                round(cation_atom_rad_std, 3), round(cation_cryst_rad_std, 3),
                round(cation_atom_rad_ptp, 3), round(cation_cryst_rad_ptp, 3),
                doc['state'], 
                sep=', ', flush=True)
        else :
            # calcuation formation energy
            total_ene_pa = doc['output']['final_energy_per_atom']
            formation_ene_pa = total_ene_pa - total_atoms_ene / (anion_num * 2.0)
            
            # check if the geometry optimzation breaks the prototype structure
            # by find if all the atoms keep the original coordinates
            # and the structure becomes somehow amphorous after relaxation
            stru = pymatgen.Structure.from_dict(doc['output']['crystal'])
            keep_stru = check_coords(structure=stru, coord_num=coord_num)
            if keep_stru:
                bond_len_mean, bond_len_std, bond_len_ptp = bond_length_statis(structure=stru)
                print(stru_type, anion, cation,  
                    round(cation_atom_rad_std, 3), round(cation_cryst_rad_std, 3),
                    round(cation_atom_rad_ptp, 3), round(cation_cryst_rad_ptp, 3),
                    doc['state'], 
                    round(formation_ene_pa, 3), 
                    round(doc['analysis']['max_force'], 3), 
                    round(doc['analysis']['delta_volume'], 3), 
                    round(doc['analysis']['percent_delta_volume'], 3),                 
                    keep_stru,                 
                    round((bond_len_mean), 3), round(bond_len_std, 3), round(bond_len_ptp, 3),
                    sep=', ', flush=True)
                if not os.path.exists(calc_dir + doc['dir_name'].split('/')[-1]):
                    os.mkdir(calc_dir + doc['dir_name'].split('/')[-1])
                    new_poscar = calc_dir + doc['dir_name'].split('/')[-1] + '/POSCAR'
                    stru.to(fmt='poscar', filename=new_poscar)
            else:
                print(stru_type, anion, cation, 
                    round(cation_atom_rad_std, 3), round(cation_cryst_rad_std, 3),
                    round(cation_atom_rad_ptp, 3), round(cation_cryst_rad_ptp, 3),
                    doc['state'], 
                    round(formation_ene_pa, 3), 
                    round(doc['analysis']['max_force'], 3), 
                    round(doc['analysis']['delta_volume'], 3), 
                    round(doc['analysis']['percent_delta_volume'], 3),                 
                    keep_stru, 
                    sep=', ', flush=True)
 