





def elements_count(data):
    '''
    Count the elements distribution histogram in the compounds dataset,  
    and the covariance matrix of each pair of elements

    Args:
        data: a list of dicts with CIFs
    Return:
        write element histogram in elem_histo file
        write element-wise covariance matrix in elem_matrix file
    '''
    # define the indice by call element_indice function
    element_indice()

    # initialize two dicts
    elem_histo = {}
    elem_matrix = {}
    for elem in periodic_table:
        elem_histo[elem] = 0
        elem_matrix[elem] = {}
        for elem2 in periodic_table:
            elem_matrix[elem][elem2] = 0

    for d in data:
        elements = Structure.from_str(d['cif'], fmt='cif').symbol_set
        for elem in elements:
            elem_histo[elem] += 1
        for elem_pair in combinations(elements, 2):
            elem_matrix[elem_pair[0]][elem_pair[1]] += 1

    for elem in periodic_table:
        for elem2 in periodic_table:
            if elem_matrix[elem][elem2] == 0 and elem != elem2:
                elem_matrix[elem][elem2] = elem_matrix[elem2][elem]

    elem_histo = {elem:count for elem,count in elem_histo.items() if count != 0}
    elem_histo = DataFrame([elem_histo])
    elem_histo = elem_histo[(elem_histo != 0)]
    elem_histo.T.to_csv('elem_histo', sep=' ')

    #
    elem_matrix = DataFrame.from_dict(elem_matrix)
    # remove columns and rows with all zero
    elem_matrix = elem_matrix.loc[:, (elem_matrix != 0).any(axis=0)] # column
    elem_matrix = elem_matrix[(elem_matrix != 0).any(axis=1)]   # row
    elem_matrix.to_csv('elem_matrix', sep=' ')

    return






def composition_similarity(baseline_id, data, index='z_number_78'):
    '''
    Calcalate the earth mover's distance between the baseline structure and all others
    in the dataset 

    Args:
        data: a pandas dataframe, index is mp-ids, and columns is element symbols
    Return:
        a pandas dataframe of pairwise distance between baseline id and all others
    '''
    # define the indice by call element_indice function
    element_indice()
    
    dist = []
    elem_similarity_file = os.path.join(sys.path[0], 'similarity_matrix.csv')
    dist_matrix = pd.read_csv(elem_similarity_file, index_col='ionA')
    dist_matrix = 1 / (np.log10(1 / dist_matrix + 1))

    if index == 'pettifor':
        dist_matrix = dist_matrix.reindex(columns=pettifor, index=pettifor) 
    # the earth mover's distance package need order C and float64 astype('float64')
    dist_matrix = dist_matrix.values.copy(order='C')

    compo_emd = pd.DataFrame(columns=[baseline_id])
    mp_id_1 = baseline_id
    for mp_id_2 in tqdm(data.index, mininterval=60):
        compo_emd.loc[mp_id_2] = emd(data.loc[mp_id_1].values.copy(order='C'), 
                                    data.loc[mp_id_2].values.copy(order='C'), 
                                    dist_matrix)
    return compo_emd


def composition_similarity_matrix(data, indice=None, index='z_number_78'):
    '''
    Calcalate pairwise earth mover's distance of all compositions, the composition should be 
    a 78-element vector, as the elemental similarity_matrix is 78x78 matrix in the order of  
    atom numbers

    Args:
        data: pandas dataframe of element vectors for all the structures
        index: see function element_indice for details
            z_number_78: (default) in the order of atomic number, this is default 
                because the similarity matrix is in this order
            z_number: see periodic_table in function element_indice
            pettifor: see function element_indice for details
            modified_pettifor: see element_indice
            elem_present: the vector only contain the elements presented in the dataset
    Return:
        a pandas dataframe of pairwise EMD with mp-ids as index
    '''
    # define the indice by call element_indice function
    element_indice()

    # if indice is None, then loop over the whole dataset
    if not indice:
        indice = [0, len(data)]
    
    dist = []
    elem_similarity_file = os.path.join(sys.path[0], 'similarity_matrix.csv')
    dist_matrix = pd.read_csv(elem_similarity_file, index_col='ionA')
    dist_matrix = 1 / (np.log10(1 / dist_matrix + 1))

    if index == 'pettifor':
        dist_matrix = dist_matrix.reindex(columns=pettifor, index=pettifor) 
    # the earth mover's distance package need order C and float64 astype('float64')
    dist_matrix = dist_matrix.values.copy(order='C')

    compo_emd = pd.DataFrame([])
    for i1 in range(indice[0], indice[1]):
        mp_id_1 = data.index[i1]
        for i2, mp_id_2 in enumerate(data.index):
            if i1 <= i2:
                emd_value = emd(data.loc[mp_id_1].values.copy(order='C'), 
                                data.loc[mp_id_2].values.copy(order='C'), 
                                dist_matrix)
                compo_emd.loc[mp_id_1, mp_id_2] = emd_value
            else:
                compo_emd.loc[mp_id_1, mp_id_2] = np.nan
    return compo_emd






def emd_of_two_compositions(y_test, y_pred, pettifor_index=True):
    '''
    Calcalate the earth mover's distance of two compositions, the composition should be 
    a 78-element vector, as the similarity_matrix is 78x78 matrix in the order of  
    atom numbers

    Args:
        y_test: test data with dimension n*78, the second dimension is by default in the  
            pettifor order, see 'composition_one_hot' function in data_explore.py
        y_pred: prediction data
        pettifor_index: whether transform the element order from the peroidic table
            number to Pettifor number
    Return:

    '''
    dist = []
    elem_similarity_file = os.path.join(sys.path[0], 'similarity_matrix.csv')
    dist_matrix = pd.read_csv(elem_similarity_file, index_col='ionA')
    dist_matrix = 1 / (np.log10(1 / dist_matrix + 1))

    if pettifor_index:
        pettifor = ['Cs', 'Rb', 'K', 'Na', 'Li', 'Ba', 'Sr', 'Ca', 'Yb', 'Eu', 'Y',  'Sc', 'Lu', 'Tm', 'Er', 'Ho', 
            'Dy', 'Tb', 'Gd', 'Sm', 'Pm', 'Nd', 'Pr', 'Ce', 'La', 'Zr', 'Hf', 'Ti', 'Nb', 'Ta', 'V',  'Mo', 
            'W',  'Cr', 'Tc', 'Re', 'Mn', 'Fe', 'Os', 'Ru', 'Co', 'Ir', 'Rh', 'Ni', 'Pt', 'Pd', 'Au', 'Ag', 
            'Cu', 'Mg', 'Hg', 'Cd', 'Zn', 'Be', 'Tl', 'In', 'Al', 'Ga', 'Pb', 'Sn', 'Ge', 'Si', 'B',  'Bi', 
            'Sb', 'As', 'P',  'Te', 'Se', 'S', 'C', 'I', 'Br', 'Cl', 'N', 'O', 'F', 'H']
        dist_matrix = dist_matrix.reindex(columns=pettifor, index=pettifor) 
    dist_matrix = dist_matrix.values

    # the earth mover's distance package need order C and float64
    dist_matrix = dist_matrix.copy(order='C')
    for y_t, y_p in zip(y_test, y_pred):
        dist.append(emd(y_t.astype('float64'), y_p.astype('float64'), 
                    dist_matrix.astype('float64')))
    
    #print(len(y_pred), np.count_nonzero(np.array(dist)))    
    return np.stack(dist)


    



