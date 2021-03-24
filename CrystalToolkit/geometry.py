
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.transformations.standard_transformations import SupercellTransformation, RotationTransformation
from pymatgen.transformations.site_transformations import TranslateSitesTransformation


def bond_length_statis(structure): 
    '''
    Calculated some bond length statistics values in the structure
    the nearest atoms pairs are determined using pymatgen CrystalNN module
    EconNN BrunnerNN_real do not work well

    Args:
        structrue: pymatgen structure
    Return:
        the mean value, standard deviation and the range of all the cation-anion
        bond length in the crystal structure
    '''
    nn = CrystalNN()
    bond_len = []
    for j in range(len(structure)) : 
        for i in nn.get_nn_info(structure,j) :       
            bond_len.append(structure[j].distance(structure[i['site_index']]))
    bond_len = np.array(bond_len)
    return bond_len.mean(), bond_len.std(), bond_len.ptp()


def change_lattice_constants(structure, lattice, scale=False):
    '''
    Change lattice constant in a pymatgen strucure

    Args:
        structure: pymatgen strucure
        lattice: a 3x3 list, like 
                    [[None,1.0,None], [None,None,None], [None,None,15]]
                    where None means this item will not be changed
        scale: is true when the lattice need to scale, i.e. latt * scale
                is False when the new lattice directly past to the structure
    Return:
        changed pymatgen structure 
    '''
    latt = structure.lattice.as_dict()
    for i in range(3):
        for j in range(3):
            if lattice[i][j] is not None:
                if scale:
                    latt['matrix'][i][j] = lattice[i][j] * latt['matrix'][i][j]
                else:
                    latt['matrix'][i][j] = lattice[i][j]

    latt = Lattice.from_dict(latt)
    structure.lattice = latt
    return structure


def check_coords(structure, coord_num):
    '''
    Find if a cell keeps the topological crystal strucure after the DFT relax
    by check the coornadinate number of all the atoms

    Args:
        structure: pymatgen structure
        coord_num: coordinate number in the ideal structure
    Return:
        A boolean that is true when the structure is kept after geometry relaxation
    '''
    nn = CrystalNN()
    keep_coord_num = True
    for atom_site in range(len(structure)) :
        if len(nn.get_nn_info(structure, atom_site)) != coord_num :
            keep_coord_num = False
            break
    return keep_coord_num


def minimum_bond_length(structure):
    '''
    Find the smallest atomic distance using pymatgen 
    structure.distance_matrix method

    Args:
        structure: pymatgen structure 
    Return:
        the value of minimum bond length
    '''
    a = structure.distance_matrix
    return a[np.nonzero(a)].min()


def multiple_structure_match(structure, structure_list):
    '''
    Check whether the input structure matches any structure in a list by
        1. first check the density differenc
        2. then use pymatgen StructureMatcher method

    Args:
        structure: pymatgen structure to compare
        structure_list: a list of structures
    Return:
        A boolean that is true when the input structre match any structure
        in the list
    '''
    is_duplic_stru = False
    if len(structure_list) == 0:
        pass
    else:
        dens1 = structure.density
        for stru_proto in structure_list:
            dens2 = stru_proto.density
            if abs(dens1 - dens2) < 0.05:
                sc = StructureMatcher()  
                if sc.fit(structure, stru_proto) :
                    is_duplic_stru = True
                    break
    return is_duplic_stru  
    
    
def planar_structure_normalization(structure):
    '''
    This function does the following:
        1. check whether the structure is planar using coordniates standard deviation
        2. move the planar layer to the center of c-direction

    Args:
        structure: pymatgen structure
    Return:
        a boolean whether the structure is planar
        tranformed pymatgen structure
    '''
    tol = 1E-3 # tolerance to check whether the structure is planar
    is_planar = True

    coords = structure.frac_coords
    ts = TransformedStructure(structure, [])
    
    if np.std(coords[:,2]) < tol : 
        center_translate = 0.5 - coords[:,2].mean()
    elif np.std(coords[:,0]) < tol :
        ts.append_transformation(SupercellTransformation([[0,0,1],[0,1,0],[1,0,0]]))
        ts.append_transformation(RotationTransformation([0,1,0], 90))        
        center_translate = 0.5 - coords[:,0].mean()
    elif np.std(coords[:,1]) < tol :
        ts.append_transformation(SupercellTransformation([[1,0,0],[0,0,1],[0,1,0]]))
        ts.append_transformation(RotationTransformation([1,0,0], 90)) 
        center_translate = 0.5 - coords[:,1].mean()
    else : 
        is_planar = False
        transformed_structure = None

    if is_planar:
        ts.append_transformation(TranslateSitesTransformation(
                            list(range(len(structure))), [0,0,center_translate]))        
        # Use pymatgen 2019.7.2, ts.structures[-1] may change in a newer version           
        transformed_structure = ts.structures[-1]
    
    return is_planar, transformed_structure


def rectangle_transform_2d(structure):
    '''
    Check whether the lattice can be transformed to a rectang lattice
    if so, return the a and b lattice constant, which will be provided
    to other function for image processing

    Args:
        structure: pymatgen structure
    Return:
        A boolean that is true when the lattice has been transformed
        a and b lattice of the rectangle
    '''
    transform_rectangle = True
    spacegrpnum = SpacegroupAnalyzer(structure, symprec=0.1).get_space_group_number()
    if spacegrpnum >= 16 and spacegrpnum <= 74 : # 16-74 is rectangler lattice
        alat = structure.lattice.matrix[0][0]
        blat = structure.lattice.matrix[1][1]
    elif spacegrpnum >= 75 and spacegrpnum <= 142 : # 75-142 is square lattice
        alat = structure.lattice.matrix[0][0]
        blat = alat
    elif spacegrpnum >= 168 and spacegrpnum <= 194 : # 168-194 is hexgonal lattice        
        alat = structure.lattice.matrix[0][0]
        blat = alat * 1.732
    else :
        transform_rectangle = False
        alat, blat = 0.0, 0.0
    return transform_rectangle, alat, blat


def symmetry_order_2d(structure, point_group_list):
    '''
    Find symmetry order in a two dimensional structure

    Args:
        structure: pymatgen structure
        point_group_list
    Return:
        Number of symmetry operations, and
        site symmetry order of the input structure
    '''
    s = SpacegroupAnalyzer(structure, symprec=0.1)
    point_group_symbols = s.get_symmetry_dataset()['site_symmetry_symbols']
    point_group_order = 0
    for point_group_symbol in point_group_symbols :
        point_group_order += point_group_list[point_group_symbol.replace('.', '')]
    # return s.get_space_group_symbol(), s.get_space_group_number(),
    return len(s.get_space_group_operations()), point_group_order/len(point_group_symbols)


             
