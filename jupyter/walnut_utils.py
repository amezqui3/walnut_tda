import tifffile as tf
from scipy import ndimage, signal
import numpy as np
import numpy.polynomial.polynomial as P
import os
import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import unionfind as UF

def normalize_density(img, npz):
    resol = 2**(img.dtype.itemsize*8)
    rescale = np.round(npz(np.arange(resol, dtype=img.dtype)))
    rescale[rescale < 0] = 0
    rescale[rescale > resol - 1] = resol - 1
    rescale.astype(img.dtype)

    with np.nditer(img, flags=['external_loop'], op_flags=['readwrite']) as it:
        for x in it:
            x[...] = rescale[x]

    return img

def get_largest_element(comp):
    labels,num = ndimage.label(comp, structure=ndimage.generate_binary_structure(comp.ndim, 1))
    print(num,'components')
    hist,bins = np.histogram(labels, bins=num, range=(1,num+1))
    argsort_hist = np.argsort(hist)[::-1]
    print(np.sort(hist)[::-1][:20])

    j = 0
    i = argsort_hist[j]
    mask = labels==i+1
    box0 = comp.copy()
    box0[~mask] = 0

    return box0

def fill_component(comp, x=True, y=True, z=True):
    rcomp = comp.copy()
    rcomp[rcomp > 0] = 1

    if x:
        for k in range(rcomp.shape[0]):
            rcomp[k,:,:] = ndimage.binary_fill_holes(rcomp[k,:,:])
        print('Closed X')
    if y:
        for k in range(rcomp.shape[1]):
            rcomp[:,k,:] = ndimage.binary_fill_holes(rcomp[:,k,:])
        print('Closed Y')
    if z:
        for k in range(rcomp.shape[2]):
            rcomp[:,:,k] = ndimage.binary_fill_holes(rcomp[:,:,k])
        print('Closed Z')

    return rcomp

def clean_zeroes(img, pad=2):
    dim = img.ndim
    orig_size = img.size

    cero = np.arange(2*dim)

    for k in range(dim):
        ceros = np.all(img == 0, axis = (k, (k+1)%dim))

        for i in range(len(ceros)):
            if(~ceros[i]):
                break
        for j in range(len(ceros)-1, 0, -1):
            if(~ceros[j]):
                break
        cero[k] = i
        cero[k+dim] = j+1
    for i in range(dim):
        cero[i] -= 2
    for i in range(dim, len(cero)):
        cero[i] += 2
    cero[cero < 0] = 0
    img = img[cero[1]:cero[4], cero[2]:cero[5], cero[0]:cero[3]]

    print(round(100-100*img.size/orig_size),'% reduction from input')

    return img, cero

def collapse_dimensions(img):
    snaps = []
    for i in range(img.ndim):
        snaps.append(np.sum(img, axis=i))
    return snaps

def collapse_dimensions_max(img):
    snaps = []
    for i in range(img.ndim):
        snaps.append(np.max(img, axis=i))
    return snaps

def plot_collapse_dimensions(snaps, bname='bname', tissue='tissue', display=False, writefig=False, dst='./'):
    fig, ax = plt.subplots(1,len(snaps),figsize=(6*len(snaps),6))
    for i in range(len(snaps)):
        ax[i].imshow(snaps[i], cmap='inferno', origin='lower');
    plt.suptitle(bname + ' ' + tissue + ' collapse', fontsize=20);
    plt.tight_layout()

    if writefig:
        filename = dst + bname + '_' + '_'.join(tissue.split(' ')) + '.jpg'
        plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight');
        if not display:
            plt.close();
