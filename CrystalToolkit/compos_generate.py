




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













