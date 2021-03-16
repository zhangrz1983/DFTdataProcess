import cv2
import numpy as np


def charge_file_crop_2d(image, structure, alat, blat, xvec=20.0, yvec=20.0, 
                        xgrid=256, ygrid=256):
    '''
    Use the lattice vectors to cut the rebundant areas
    i.e. only the periodic lattice is kept, and put in the center of the image
    Using white color padding other areas to make a 256x256 resolution
    
    Args:
        alat, blat：lattice constants of the 2D cell
        xvec, yvec: the 2D lattice vectors of the crystal cell
        xgrid, ygrid: pixel resolution of the charge density image
    Return:
        A cropped image with white padding near the edges
    '''
    WHITE = [0, 0, 0] # color code used for padding
    xreso =  xvec / xgrid  # charge density resolution when doing slice
    yreso =  yvec / ygrid  
    na = int(xvec / alat)
    nb = int(yvec / blat)
    aposit = int(alat * na / xreso)
    bposit = int(blat * nb / yreso)
    if ( aposit % 2 ) != 0 :
        aposit = aposit + 1
    if ( bposit % 2 ) != 0 :
        bposit = bposit + 1
    apad = int(( xgrid - aposit ) / 2 )
    bpad = int(( ygrid - bposit ) / 2 )

    # mind that height goes first in the crop function
    image_crop = image[0:bposit, 0:aposit] 
    new_image = cv2.copyMakeBorder(image_crop, bpad, bpad, apad, apad,  
                        cv2.BORDER_CONSTANT, value=WHITE)
    return new_image


def find_latt_vectors_contour_2d(image, latt_thresh_value=1, round_value=1):
    '''
    Find the two lattice vectors using OpenCV contour method 
    the threshold in the cv2.threshold is very low, to give the whole area of the lattice
    the lattice area should the largest contour found by cv2.RETR_TREE
    
    Alternatively the following code can be used, but RETR_EXTERNAL might not work well
        cnts_latt = cv2.findContours(thresh_latt.copy(), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)[0]
        return cv2.boundingRect(cnts_latt[0])

    Args:
        latt_thresh_value: determine which region is turned into black in a binary image
                            usually this value is set to be low (e.g. 1)
        round_value : the OpenCV detectd areas is somehow larger than the real lattice vectors, 
                        possible reasons include:
                            when put the cropped image in the center of 256*256 resolution,
                            structures with odd pixel number are rounded to even number by adding 1    
    Return:
        x_start, y_start, x_latt, y_latt : which define the rectangle of the lattice vectors
    '''
    thresh_latt = cv2.threshold(img, latt_thresh_value, 255, cv2.THRESH_BINARY)[1]
    cnts_latt = cv2.findContours(thresh_latt.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[0]
    areas = [cv2.contourArea(c) for c in cnts_latt]
    max_index = np.argmax(areas)
    max_cnt = cnts_latt[max_index]
    x_start, y_start, x_latt, y_latt = cv2.boundingRect(max_cnt)
    x_start = x_start + round_value
    y_start = y_start + round_value
    x_latt = x_latt - round_value * 2
    y_latt = y_latt - round_value * 2
    return x_start, y_start, x_latt, y_latt


def find_atom_coords_contour_2d(img, x_start, y_start, x_latt, y_latt,
                                atom_thresh_value=63, 
                                max_atom_rad=10.0, min_atom_rad=3.0):
    '''
    find the 2D coordinates of the atoms using OpenCV contour method
    
    Args:
        atom_thresh_value: determine which region is turned into white in a binary image
            usually this value is set to be low (i.e. 15) to aviod the region with low 
            electronic density to be regarded as atom. Note that Binary_inv is used for 
            cv2.threshold, so void area (no electrons) are WHITE and the blank area near 
            the edge is filled with BLACK to detect atoms at the edges to distwish atoms 
            and void area, use the raduis of the circles
        min_atom_rad, max_atom_rad: criterion for atom detection using atomic raduis
    Returns:

    '''
    stencil = np.zeros(img.shape).astype(img.dtype)
    stencil = cv2.rectangle(stencil,
                            (x_start, y_start), (x_start + x_latt, y_start + y_latt),
                            (255, 255, 255), -1 )
    thresh_atom = cv2.threshold(img, atom_thresh_value, 255, cv2.THRESH_BINARY_INV)[1]
    # use a stencil to add white color to the void area
    thresh_atom = cv2.bitwise_and(thresh_atom, stencil)
    cnts_atom = cv2.findContours(thresh_atom.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[0]

    coords = []
    for cnt in cnts_atom:
        if ( cv2.minEnclosingCircle(cnt)[1] < max_atom_rad and
                cv2.minEnclosingCircle(cnt)[1] > min_atom_rad ):
            coords.append(cv2.minEnclosingCircle(cnt)[0])
    return coords


def merge_period_cond_atoms_2d(img, coords, x_start, y_start, x_latt, y_latt, 
                            symm_crit=5, edge_detect=6.0):
    '''
    Merge the atoms due to periodic conditions, first x then y
    during image generation there might be redundant region, remove them using x/y_pero_cond
    so x_latt and y_latt may change

    Args:
        edge_detect : if the distance between the atom center and the edge is smaller than
            edge_detect, then this atom is considered to be on the edge and subject to 
            the periodic condition
        symm_crit : if the pixel distance between two atoms (usually at image edges) is
            smaller than symm_crit, the two atoms are considered as one under periodic
            condition
        x_pero_cond, y_pero_cond : if the two atoms near the edge are detemined to be one 
            atom, then the distance between them (x_pero_cond*2 or y_pero_cond*2)
            should be set to 0 to make them one atom, and the lattice vector should minus 
            this value
    '''
    x_pero_cond, y_pero_cond = 0.0, 0.0
    for coord in coords[:] :
        if coord[0] - x_start < edge_detect : 
            for coord2 in coords :
                # whether the two atoms are reflective across the x center line and have same y coordinate
                if (abs( (coord[0] + coord2[0]) / 2 - img.shape[0] / 2 ) < symm_crit  and
                    abs(coord[1] - coord2[1]) < symm_crit ):
                    x_pero_cond = (coord[0] - x_start) * 2
                    coords.remove(coord)
    for coord in coords[:] :
        if coord[1] - y_start < edge_detect : 
            for coord2 in coords :
                if (abs( (coord[1] + coord2[1]) / 2 - img.shape[1] / 2 ) < symm_crit  and
                    abs(coord[0] - coord2[0]) < symm_crit ):
                    y_pero_cond = (coord[1] - y_start) * 2
                    coords.remove(coord)
    x_latt = x_latt - x_pero_cond
    y_latt = y_latt - y_pero_cond
    return x_latt, y_latt, coords


def write_image_to_poscar_2d(image_file, x_start, y_start, x_latt, y_latt, coords,
                            xvec=20.0, yvec=20.0, xgrid=256, ygrid=256):
    '''
    Args:
    Return:
        pymatgen structure object
    '''
    with open (image_file.replace('png', 'vasp'), 'w+') as f :
        f.write(image_file + '\n')
        f.write('1.0 \n')
        f.write(str(x_latt / xgrid * xvec) + ' 0.0  0.0 \n')    
        f.write('0.0  ' + str(y_latt / ygrid * yvec) + '  0.0 \n')
        f.write('0.0 0.0 15.0 \n')
        f.write(image_file.split('_')[0] + ' \n') # print the element name
        f.write(str(len(coords)) + '\n')
        f.write('D \n')

        for coord in coords :
            f.write(str( (coord[0] - x_start) / x_latt ) + '  ' + 
                    str( (coord[1] - y_start) / y_latt ) + ' 0.5 \n')


