
import os
import glob
import numpy as np
import cv2

import image_file_process

if __name__ = '__main__':
    for filein in glob.glob('*.png'):
        img = cv2.imread(filein)
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        x_start, y_start, x_latt, y_latt = image_file_process.find_latt_vectors_contour_2d(img=gray)
        coords = image_file_process.find_atom_coords_contour_2d(img=gray,
                                                                x_start=x_start, y_start=y_start,
                                                                x_latt=x_latt, y_latt=y_latt)
        x_latt, y_latt, coords = image_file_process.merge_period_cond_atoms_2d(img=gray, coords=coords, 
                                                                                x_start=x_start, y_start=y_start,
                                                                                x_latt=x_latt, y_latt=y_latt)
        image_file_process.write_image_to_poscar_2d(image_file=filein,
                                                    x_start=x_start, y_start=y_start,
                                                    x_latt=x_latt, y_latt=y_latt, coords=coords)

    




