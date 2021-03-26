






def bonding_type(structure):
    '''
    Calculate bonding in a given structure.

    Args:
        structure: a pymatgen structure
    Return:
        A list of atomic pairs which form bonding
    '''
    nn = CrystalNN()
    bond_elem_list = []
    bond_num_list = []
    for i in list(range(len(struct))):
        site1 = struct[i].species_string
        num1 = struct[i].specie.number
        for neigh in nn.get_nn_info(struct, i):
            bond_elem_list.append(' '.join(sorted([site1, neigh['site'].species_string])))
            bond_num_list.append(' '.join(list(map(str, sorted([num1, neigh['site'].specie.number])))))

    bond_elem_list = list(set(bond_elem_list))
    bond_num_list = list(set(bond_num_list))
   
    return bond_elem_list, bond_num_list
            

def bonding_matrix(data):
    '''
    Convert the atomic pair of bonding into a matrix representation
    Remove the elements that don't exsit, then flatten it.

    Args:
        data: a list of dicts with CIFs
    Return:
        A 2D vector, first dimension is number of samples.
    '''
    # define the indice by call element_indice function
    element_indice()
    
    all_bond_matrix = []
    num_elements = len(periodic_table_78)
    zeor_bond_matrix = pd.DataFrame(np.zeros((num_elements, num_elements)), 
                                    index=periodic_table_78, columns=periodic_table_78)
    for d in data:
        bond_matrix = zeor_bond_matrix
        bond_list = d['bond_elem_list']
        for bond in bond_list:
            elem1, elem2 = bond.split()
            bond_matrix.loc[elem1, elem2] = 1
            bond_matrix.loc[elem2, elem1] = 1
        all_bond_matrix.append(bond_matrix.values.flatten())
    # HERE NEED TO DELETE ZERO ROW/COLUMN AND THEN FLATTEN
    return np.stack(all_bond_matrix)




def average_coordination(structure):
    '''
    Calculation of average coordination number over every site in a structure
    using Vorini method

    Args:
        structure: pymatgen structure
    Return:
        Average coordination number
    '''
    nn = CrystalNN()
    ave_coord_num = []
    for atom_site in range(len(structure)):
        ave_coord_num.append(len(nn.get_nn_info(structure, atom_site)))
    return np.array(ave_coord_num).mean()


def bond_stat_per_site(structure):
    '''
    Get the bond length standard deviation of each site in the unit cell, the average
    Get the coordination number of each site, then standard deviation

    Args:
        structure: pymatgen structure
    Return:
        Average of bond length std of each site
        Standard deviation of coordination number
    '''
    nn = CrystalNN()
    coord_num = []
    bond_len_std = []
    for atom_site in range(len(structure)):
        struct_nn = nn.get_nn_info(structure, atom_site)
        coord_num.append(len(struct_nn))
        bond_len = []
        for i in struct_nn:       
            bond_len.append(structure[atom_site].distance(structure[i['site_index']]))
        bond_len_std.append(np.array(bond_len).std())
    return round(np.array(bond_len_std).mean(), 3), \
            round(np.array(coord_num).std(), 3)



def bond_to_atom(data, nelem=78):
    '''
    Convert atomic pairs bonding into the exsitance of elements

    Args:
        Data: y_data of flatten matrix
        nelem: number of elements in the periodic table, determines the shape
            of the bonding matrix
    Return:
        One-hot vector of elements
    '''
    new_data = []
    for y in data:
        y = np.reshape(y, (nelem, nelem))
        # make all the non-zero values to 1
        # NOTETHAT the sum method is only valid for a symmetric bonding matrix
        new_data.append(np.sum(y, axis=0).astype(bool).astype(int))
    return np.stack(new_data)









