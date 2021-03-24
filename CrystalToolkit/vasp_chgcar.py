
import sys
import os
import linecache
import numpy as np
from PIL import Image


def chgcar_file_slice_2d(density, atoms, xvec=20.0, yvec=20.0, xgrid=256, ygrid=256):
    '''
    Modified from 
        slice.py - version 0.3 
        Copyright (C) 2008, 2009 Matthew Dyer and Jonas Bjork
    Main change: 
        1. Python 2 to python 3, and remove the interactive mode. 
        2. Add lattice start point to cover the whole 2D lattice
        3. Change the output format to matrix format
    
    Args:
        density: ASE charge density object 
        atoms: ASE crystal structure object
    Return:
        nplotgrid: number of points of the 2D charge density
        xarray, yarray: 
        density2D: 2D change density martrix
    '''
    xvector = np.array([xvec, 0.0, 0.0])
    yvector = np.array([0.0, yvec, 0.0])
    nplotgrid = np.array([xgrid, ygrid])

    # find lattic start point, to include the whole lattice in the sliced density file
    alatx, alaty = linecache.getline(filein, 3).split()[0:2]
    blatx, blaty = linecache.getline(filein, 4).split()[0:2]
    clat = float(linecache.getline(filein, 5).split()[2])
    xstart = min(float(alatx), float(blatx), 0.0)
    ystart = min(float(alaty), float(blaty), 0.0)
    origin = np.array([xstart, ystart, 0.5*clat])

    ngridpts = np.array(density.shape)
    totgridpts = ngridpts.prod()
    unit_cell = atoms.get_cell()

    fdensity = np.fft.fftn(density)
    cartgrid = np.zeros((nplotgrid[0], nplotgrid[1], 3), np.float)
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            cartgrid[i][j] = origin \
                             + float(i)/(nplotgrid[0] - 1) * xvector \
                             + float(j)/(nplotgrid[1] - 1) * yvector

    scaledgrid = []
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            scaledgrid.append(np.dot(cartgrid[i][j], np.linalg.inv(unit_cell))%1)
    scaledgrid = np.array(scaledgrid)
    scaledgrid = 2.*np.pi*1.j*scaledgrid

    gvectors = []
    temp = []
    for i in range(int(-ngridpts[0]/2) + 1, int(ngridpts[0]/2) + 1):
        for j in range(int(-ngridpts[1]/2) + 1, int(ngridpts[1]/2) + 1):
            for k in range(int(-ngridpts[2]/2) + 1, int(ngridpts[2]/2) + 1):
                gvectors.append([float(i), float(j), float(k)])
                temp.append(fdensity[i][j][k])
    gvectors = np.array(gvectors)
    gvectors = gvectors.T
    fdensity = np.array(temp)
    del temp

    density2D = np.zeros((nplotgrid[0] * nplotgrid[1]), np.float)
    for i in range(nplotgrid[0] * nplotgrid[1]):
        temp = np.dot(scaledgrid[i], gvectors)
        temp = np.exp(temp)
        density2D[i] = np.dot(fdensity, temp.T).real
    del temp
    del gvectors
    del fdensity

    density2D = density2D / totgridpts
    density2D = density2D.reshape(nplotgrid[1], nplotgrid[0])
    xlength = np.sqrt(np.dot(xvector, xvector.T))
    ylength = np.sqrt(np.dot(yvector, yvector.T))
    yontox = np.dot(xvector, yvector.T)/xlength
    ynormal = np.cross(xvector, yvector.T)/xlength
    ynormal = np.sqrt(np.dot(ynormal, ynormal.T))

    xarray = np.zeros((nplotgrid[1], nplotgrid[0]), np.float)
    yarray = np.zeros((nplotgrid[1], nplotgrid[0]), np.float)
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            xarray[j][i] = float(i) / float(nplotgrid[0] - 1) * xlength \
                           + float(j) / float(nplotgrid[1] - 1) * yontox
            yarray[j][i] = float(j) / float(nplotgrid[1] - 1) * ynormal

    # some values are negative in the Frontier transform in 
    # so make them 0 to avoid error when change the number matrix to image
    density2D = density2D.clip(min=0)
    
    return nplotgrid, xarray, yarray, density2D


def charge_file_2d_reformat(filein, xgrid=256, ygrid=256):
    '''
    To be deprecated
    This funciton is only for the orginal slice.py, 
    where the chgcar_file_slice_2d output three columns, i.e. x_posit, y_posit, value, 
    and these are seperated by empty lines
    This function reformat it, only keep values in the MATRIX format
    The loadtxt function automatically remove empty lines
    '''
    # failed to use numpy arrary(npy format) in the PyTorch-GAN-zoo
    # so use png format as input to PyTorch 
    #file_npy = filein.replace('chgtile', 'npy')
    file_png = filein.replace('chgtile', 'png')

    chgden_2d = np.loadtxt(filein, usecols=(2), unpack=True)
    chgden_2d = np.resize(chgden_2d, (xgrid, ygrid))
    # some values are negative in the Frontier transform in 
    # so make them 0 to avoid error when change the number matrix to image
    chgden_2d = chgden_2d.clip(min=0)
    # type must be unit 8 otherwise PIL will complain
    chgden_2d = ( 255.0 / chgden_2d.max() * 
                    (chgden_2d - chgden_2d.min()) ).astype(np.uint8)
    #np.save(file_npy, chgden_2d)
    Image.fromarray(chgden_2d).save(file_png)