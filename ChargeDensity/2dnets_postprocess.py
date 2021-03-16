'''
usage:
python 2dnets_postprocess.py CHGCAR_FILE_NAME
VASP chgcar file is required
'''

import sys
import json
import numpy as np
from pymatgen.io.vasp.outputs import Chgcar
from ase.calculators.vasp import VaspChargeDensity
from PIL import Image

from ChargeDensity.image_file_process import charge_file_crop_2d
from CrystalToolkit.geometry_properties import symmetry_order_2d, rectangle_transform_2d
from CrystalToolkit.vasp_chgcar import chgcar_file_slice_2d, charge_file_2d_reformat


file_chgcar = sys.argv[1]
stru = Chgcar.from_file(file_chgcar).poscar.structure
'''
# find the symmetry order of the 2D structure
with open('data.json') as json_file:
    symmetry_data = json.load(json_file)
point_group_list = symmetry_data['point_group_list']
symmetry_order_2d(structure=stru, point_group_list=point_group_list)
'''
# slice the vasp chgcar file into 2D
vasp_charge = VaspChargeDensity(filename=file_chgcar)
density = vasp_charge.chg[-1]
atoms = vasp_charge.atoms[-1]
*_, chgden_2d = chgcar_file_slice_2d(density=density, atoms=atoms)
# change to RGB values, type must be unit 8 otherwise PIL will complain
chgden_2d = ( 255.0 / chgden_2d.max() * 
                (chgden_2d - chgden_2d.min()) ).astype(np.uint8)
file_png = file_chgcar.replace('chgcar', 'png')
Image.fromarray(chgden_2d).save(file_png)
#cv2.imwrite(file_png, chgden_2d)

#charge_file_2d_reformat(filein=file_chgcar.replace('chgcar', 'chgtile'))

if_transform, alat, blat = rectangle_transform_2d(structure=stru)
if if_transform:
    im = cv2.imread(file_png, 0)
    new_im = charge_file_crop_2d(image=im, structure=stru, alat=alat, blat=blat)
    fileout = file_chgcar.replace('chgcar', 'padding.png')
    cv2.imwrite(fileout, new_im)

