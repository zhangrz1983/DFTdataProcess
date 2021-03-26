



def element_indice():
    '''
    type of element indice used 
    '''
    global modified_pettifor, pettifor, periodic_table, periodic_table_78

    # 78 elements without nobal gas and rare earch Ac
    # Put the elements with similar propeties together to give a better representation, see
    # Pettifor, D. G. (1990). "Structure maps in alloy design." 
    # Journal of the Chemical Society, Faraday Transactions 86(8): 1209.
    pettifor = [
        'Cs', 'Rb', 'K', 'Na', 'Li', 
        'Ba', 'Sr', 'Ca', 
        'Yb', 'Eu', 'Y',  'Sc', 'Lu', 'Tm', 'Er', 'Ho', 
        'Dy', 'Tb', 'Gd', 'Sm', 'Pm', 'Nd', 'Pr', 'Ce', 'La', 
        'Zr', 'Hf', 'Ti', 'Nb', 'Ta', 'V',  'Mo', 'W',  'Cr', 'Tc', 'Re', 
        'Mn', 'Fe', 'Os', 'Ru', 'Co', 'Ir', 'Rh', 'Ni', 'Pt', 'Pd', 'Au', 'Ag', 'Cu', 
        'Mg', 'Hg', 'Cd', 'Zn', 'Be', 'Tl', 'In', 'Al', 'Ga', 'Pb', 'Sn', 'Ge', 'Si', 'B',  
        'Bi', 'Sb', 'As', 'P',  'Te', 'Se', 'S', 'C', 'I', 'Br', 'Cl', 'N', 'O', 'F', 'H'
    ]

    # 103 elements, A modified version of Pettifor, see
    # Glawe, H., et al. (2016). New Journal of Physics 18(9): 093011.
    # "The optimal one dimensional periodic table: a modified Pettifor chemical scale from data mining." 
    modified_pettifor = [
        'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 
        'Fr', 'Cs', 'Rb', 'K', 'Na', 'Li', 'Ra', 'Ba', 'Sr', 'Ca', 
        'Eu', 'Yb', 'Lu', 'Tm', 'Y', 'Er', 'Ho', 'Dy', 'Tb', 'Gd', 'Sm', 'Pm', 'Nd', 'Pr', 'Ce', 'La', 
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 
        'Sc', 'Zr', 'Hf', 'Ti', 'Ta', 'Nb', 'V', 'Cr', 'Mo', 'W', 'Re', 
        'Tc', 'Os', 'Ru', 'Ir', 'Rh', 'Pt', 'Pd', 'Au', 'Ag', 'Cu', 
        'Ni', 'Co', 'Fe', 'Mn', 'Mg', 'Zn', 'Cd', 'Hg', 
        'Be', 'Al', 'Ga', 'In', 'Tl', 'Pb', 'Sn', 'Ge', 'Si', 'B', 'C', 
        'N', 'P', 'As', 'Sb', 'Bi', 'Po', 'Te', 'Se', 'S', 'O', 'At', 'I', 'Br', 'Cl', 'F', 'H'
    ]

    # 89 elements in the order of atomic Z number
    # This is the periodic table which DFT pesudopotentials are avaiable
    periodic_table = [
        'H',  'He', 
        'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 
        'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 
        'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
        'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe', 
        'Cs', 'Ba', 
                    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
                          'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 
    ]

    # 78 elements version by removing nobal gas and six rare-earth Ac-row elements
    periodic_table_78 = [
        'H',  
        'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  
        'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl',  
        'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 
        'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  
        'Cs', 'Ba', 
                    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
                          'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 
    ]







def similarity_matrix(input_file='dist_matrix', normalize='inverse', order='pt_number'):
    '''
    Read the original file from the below reference, and convert into a dictionary

    Hautier, G., et al. (2011). 
    "Data mined ionic substitutions for the discovery of new compounds." 
    Inorganic Chemistry 50(2): 656-663.

    Args:
        input_file: the preprocessed supporting txt file (from the above paper)
                using the shell script below
        normalize: method for value normalization
            bound: all value divide by 20 (based on the maximum value)
            log: all value got log10
            inverse: 1 / values
        order: the order of element in the matrix
            pt_number: the order in periodic table, i.e. atomic number
                    typically for calculation purpose
            pettifor: typically for visualization purpose, 
     Return:
        a pandas dataframe of the similarity matrix
    ====================================================
    #!/bin/sh
    outfile=dist_matrix
    > $outfile
    for i in `seq 78` # up to Bi and remove nobel gas
    do
        e1=`sed -n "$i p" mendlev` # element symbols in Pettifor order
        p1=`sed -n "$i p" Pettifor` # element symbol + valency in Pettifor order
        for j in `seq 78`
        do
            e2=`sed -n "$j p" mendlev`
            p2=`sed -n "$j p" Pettifor`
            if [ $i -gt $j ]
            then
                r=`grep $p1 ic102031h_si_001.txt | grep $p2`
                if [ -z "$r" ]
                then
                    grep -w $e1 ic102031h_si_001.txt | grep -w $e2 | head -n 1 >> $outfile
                else
                    echo $r >> $outfile
                fi
            fi
        done
    done
    sed -i 's/:/ /g' $outfile # make the valency a seperate column
    =========================================================
    '''
    # define the indice by call element_indice function
    element_indice()

    # note that the index_col is after the use of selection
    d = pd.read_csv(input_file, sep=' ', index_col=[0,1], usecols=[0,2,4])
    # make it a two dimensional matrix
    d = d.unstack()
    # drop the multilevel when ustack
    d.columns = d.columns.droplevel()

    if order == 'pt_number':
        index = periodic_table_78
    elif order == 'pettifor':
        index = pettifor
    # reindex, i.e. change the order of column and rows to Pettifor order  
    d = d.fillna(0)
    d = d.reindex(columns=index, index=index, fill_value=0) 
    d = d + d.transpose()
    # the maximum similarity number is 18.6, set the same element to 20, a bit arbitary
    np.fill_diagonal(d.values, 20)
    # then fill the zeros
    d.loc[:,'Pm'] = d.loc[:,'Sm'] # column, same as:  d['Pm'] = d['Sm']
    d.loc['Pm',:] = d.loc['Sm',:] # row
    # other zero mean very dissimilar, so set to a very small value
    d.replace(0, 0.1, inplace=True)

    if normalize == 'bound':
        d = d / 20
    elif normalize == 'log':
        d = np.log10(d)
    elif normalize == 'inverse':
        d = 1 / d
        np.fill_diagonal(d.values, 0)
    else:
        print('normalization method not supported')

    return d





def elements_selection(data, elem_list, mode='include'):
    '''
    Select a subset contains or not contain certain elements.

    Args:
        data: a list of dicts with CIFs
        elem_list: a list of the elements of interest or no-interest
        mode:
            include: select structures have elements in elem_list
            exclude: drop structures have elements in elem_list
            consist: select structures made up of elements in elem_list
    Return:
        A new dataset after selection
    '''
    for d in data[:]:
        elements = Structure.from_str(d['cif'], fmt='cif').symbol_set
        if mode == 'include':
            if set(elem_list).isdisjoint(elements):
                data.remove(d)
        elif mode == 'exclude':
            if not set(elem_list).isdisjoint(elements):
                data.remove(d)
        elif mode == 'consist':
            if not set(elements).issubset(elem_list):
                data.remove(d)
    return data









