

def composition_one_hot(data, method='percentage', index='z_number_78', only_elem_present=False):
    '''
    Make the composition to fix size vector like one hot array
    Note that this function is ONLY suitable for formula with integer number
    of each atom, i.e. a fomula with fractional occuption is not supported

    !!!REINDEX FOR ONLY_TYPE=FASLE NOT IMPLEMENT!!!

    Args:
        data: json dataset of structures from Materials Project
        Method: how the number of atoms is encoded
            percentage: (default) the number of atoms given as percentage, 
                i.e. the sum of the whole formula is 1
            formula: numbers in formula
            only_type: the value of each atom is set to 1
        index: see function element_indice for details
            z_number_78: (default) in the order of atomic number, this is default 
                because the similarity matrix is in this order
            z_number: see periodic_table in function element_indice
            pettifor: see function element_indice for details
            modified_pettifor: see element_indice
            elem_present: the vector only contain the elements presented in the dataset
        only_element_present: Used when the vector is reindexed. If true, only keep the 
            elements presented in the dataset
    Return:
        a pandas dataframe, index is mp-ids, and columns is element symbols
        elem_symbol is a list of elements present in the dataset, in Z number order
    '''
    # define the indice by call element_indice function
    element_indice()

    pero_tab_nums = []
    mp_index = []
    for d in data:
        struct = Structure.from_str(d['cif'], fmt='cif')
        # the struct.species method give a list of each site, 
        # e.g. for SrTiO3 the output is 
        # [Element Sr, Element Ti, Element O, Element O, Element O]
        pero_tab_nums.append([ x.number for x in struct.species ])
        mp_index.append(d['task_id'])
    
    # use the counter method to one-hot the element numbers
    # considering the output of species method above
    # the data_np is typically using " method == 'formula' "
    elem_vectors = pd.DataFrame([Counter(x) for x in pero_tab_nums])
    elem_vectors = elem_vectors.fillna(0).sort_index(axis=1)

    if method == 'percentage':
        # divide by total number of atoms to the sum will be 1
        # sum all columns in each row, divide row-wise
        elem_vectors = elem_vectors.div(elem_vectors.sum(axis=1), axis=0)
    elif method == 'only_type':
        # set all the non-zero values to 1
        elem_vectors[elem_vectors != 0] = 1

    # get the order of elements in the list and element symbols
    elem_numbers = elem_vectors.columns.tolist()
    # this gives all element symbols present in the dataset
    elem_symbols = [ Element.from_Z(x).name for x in elem_numbers ]
    
    elem_vectors.columns = elem_symbols
    elem_vectors.index = mp_index

    # Now the vectors in data_np is in the order of Z number but only
    # with the elements presented in the dataset
    # we may want to reindex the data_np
    # Note the index here is accutely column names in pandas, not pandas index
    if index != 'elem_present': 
        if index == 'z_number_78':
            elem_vectors = elem_vectors.reindex(columns=periodic_table_78)
        elif index == 'z_number':
            elem_vectors = elem_vectors.reindex(columns=periodic_table)
        elif index == 'pettifor':
            elem_vectors = elem_vectors.reindex(columns=pettifor)
        elif index == 'modified_pettifor':
            elem_vectors = elem_vectors.reindex(columns=modified_pettifor)
            
        if only_elem_present:
            elem_vectors = elem_vectors.dropna()
        else:
            elem_vectors = elem_vectors.fillna(0)

    return elem_vectors, elem_symbols





