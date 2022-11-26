import tifffile as tf
from scipy import ndimage
from scipy import spatial
from scipy import signal

import numpy as np
import pandas as pd

import os
import glob
import argparse

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import walnut_utils as wnut

thr = 1
fs = 15
resol = 0.0759
NNN = 20

mborder = ndimage.generate_binary_structure(3, 1).astype(int)
mborder[mborder == 1] = -1
mborder[1,1,1] = -np.sum(mborder) - 1

parser = argparse.ArgumentParser(description='Normalize density values of walnuts')

parser.add_argument('src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')
parser.add_argument('dst', metavar='clean_img_dst', type=str, help='path to store clean images')
parser.add_argument('bname', metavar='scan_id', type=str, help='walnut batch scan id')

args = parser.parse_args()

# src = '../clean2/'
# dst = '../rotated/'
# bname = '2014SBa_R1_T25'

dst = args.dst
src = args.src
bname = args.bname

walnut_files = sorted(glob.glob(src + bname + '/*.tif'))

wdst = dst + bname + '/'
if not os.path.isdir(wdst):
    os.makedirs(wdst)
print(wdst)

for widx in range(len(walnut_files)):

    pa, fi = os.path.split(walnut_files[widx])
    fname = os.path.splitext(fi)[0]
    print(fname)

    img = tf.imread(walnut_files[widx])
    mxval = 2**(8*img.itemsize)

    nutvol = np.sum(img > 0) * (resol ** 3)

    bottom = img[:NNN,:, :]
    top = img[-NNN:, :, :]

    topsum = np.sum(top > 0)
    botsum = np.sum(bottom > 0)

    if botsum > topsum:
        coords = np.asarray(np.nonzero(top))
        tipvox = wnut.find_tip_max(coords, 2,1,0)
        tipvox0 = tipvox.copy()
        tipvox[0] += (img.shape[0] - NNN)
    else:
        coords = np.asarray(np.nonzero(bottom))
        tipvox = wnut.find_tip_min(coords, 2,1,0)
        tipvox0 = tipvox.copy()

    bimg = img.copy().astype(int)
    bimg[bimg > 0]  = 1
    border = ndimage.convolve(bimg, mborder, mode='constant', cval=0)
    border[border < 0] = 0
    nutarea = np.sum(border) * (resol ** 2)
    border[border > 0] = 1
    nutvoxarea = np.sum(border)

    datapoints = np.asarray(np.nonzero(border))
    wmean = np.mean(datapoints, axis = 1)

    datapoints = datapoints - wmean.reshape(-1,1)
    datapoints *= resol

    fig,ax = wnut.plot_3Dprojections(datapoints, fname + ' original', alpha=0.01, writefig=True, dst=wdst, display=False);

    tipvox = tipvox - wmean
    tipvox *= resol

    flipx = np.eye(3)
    rotX, rotY, rotZ = False, False, False

    if tipvox[0] < 0:
        rotX = True
        datapoints[0] *= -1.
        tipvox[0] *= -1.
        flipx[0,0] *= -1

    thetay = np.sign(tipvox[1]) * np.arccos(tipvox[0]/np.sqrt(tipvox[0]**2 + tipvox[2]**2))
    thetaz = np.sign(tipvox[2]) * np.arccos(tipvox[0]/np.sqrt(tipvox[0]**2 + tipvox[1]**2))

    roty = np.array([[np.cos(thetay), 0,  np.sin(thetay)],
                     [0, 1, 0],
                     [-np.sin(thetay), 0, np.cos(thetay)]])

    rotz = np.array([[np.cos(thetaz), -np.sin(thetaz), 0],
                     [np.sin(thetaz),  np.cos(thetaz), 0],
                     [0, 0, 1]])

    rots = [rotz.T @ roty, rotz @ roty.T, rotz.T @ roty, rotz @ roty]
    foo = np.zeros(len(rots))

    for i in range(len(rots)):
        bar = rots[i] @ tipvox
        foo[i] = np.sum(bar[1:] ** 2)

    rotyz = rots[np.argmin(foo)]

    rcoords = rotyz @ datapoints
    rtipvox = rotyz @ tipvox

    acoords = np.abs(rcoords)
    tcoords = np.flip(rcoords[1:, (acoords[0] < thr)])

    V = wnut.ell_algebraic_fit_2d(*tcoords)
    pdict = wnut.get_ell2d_params_from_vector(V)

    Q = pdict['rot'].T @ (tcoords - pdict['origin'].reshape(-1,1))
    distances = np.zeros((Q.shape[1], 4))

    for i in range(Q.shape[1]):
        distances[i] = wnut.d_ell_point(Q[0,i], Q[1,i], pdict)

    P = pdict['rot'] @ (distances[:, 2:].T + pdict['origin'].reshape(-1,1))

    phi = np.arccos(distances[:,2]/pdict['axes'][0])
    phi[distances[:,3] < 0] = 2*np.pi - phi[distances[:,3] < 0]

    stack = wnut.even_space_ell0(721, *pdict['axes'])
    idxs = np.digitize(phi, bins=stack[0])

    sumdist = np.zeros(len(stack[0])-1)
    for i in range(len(sumdist)):
        sumdist[i] = np.mean(distances[idxs == i+1, 0]*distances[idxs == i+1, 1])
    sumdist[sumdist < 0] = 0

    cval = np.max(np.abs(tcoords))/np.max(sumdist)

    sumdist3 = np.hstack((sumdist,sumdist,sumdist))
    rawpeakidx, _ = signal.find_peaks(sumdist3, distance=300)
    foo = rawpeakidx[(rawpeakidx > len(sumdist)) & (rawpeakidx < 2*len(sumdist))]
    srtpeakidx = foo[np.argsort(sumdist3[foo])[-2:][::-1]] - len(sumdist)

    thetax = np.mean(stack[0][srtpeakidx]) - np.pi*.5 + pdict['theta']

    if np.argmax(sumdist) < len(sumdist)//2:
        thetax += np.pi
        print('rot')

    rotx = np.array([[1,0,0],
                     [0, np.cos(thetax), -np.sin(thetax)],
                     [0, np.sin(thetax),  np.cos(thetax)]])

    rotxyz = rotx @ rotyz
    bulgerot2d = rotx[1:,1:]

    fig, ax = plt.subplots(1,2,figsize=(12,6), sharex=True, sharey=True)

    ax[0].scatter(*tcoords, s=1)#, c=distances[:,0]*distances[:,1], cmap='plasma')
    ax[0].scatter(*P, s=15, c=distances[:,0]*distances[:,1], cmap='plasma')
    ax[0].plot(cval*sumdist*np.cos(stack[0][:-1] + pdict['theta']), cval*sumdist*np.sin(stack[0][:-1] + pdict['theta']), c='r')

    ax[1].scatter(*(bulgerot2d.T @ tcoords), s=1)#, c=distances[:,0]*distances[:,1], cmap='plasma')
    ax[1].scatter(*(bulgerot2d.T @ P), s=15, c=distances[:,0]*distances[:,1], cmap='plasma')

    for i in range(len(ax)):
        ax[i].axvline(0, c= 'lightgray')
        ax[i].axhline(0, c= 'lightgray')
        ax[i].set_aspect('equal');

    fig.suptitle(fname + ' Seal slice alignment', fontsize=20)
    fig.tight_layout()

    filename = wdst + fname + '_seal.jpg'
    plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

    rcoords = rotxyz @ datapoints
    rtipvox = rotxyz @ tipvox

    minis = np.min(rcoords, axis=1)
    maxes = np.max(rcoords, axis=1)

    fig, ax = wnut.plot_3Dprojections(rcoords, fname + '+rotXYZ', alpha=0.01, writefig=True, dst=wdst, display=False);

    hull = spatial.ConvexHull(rcoords.T)
    align = np.vstack((wmean, rotx, roty, rotz, rotxyz @ flipx,
                       np.array([rotX,rotY,rotZ]).astype(int),
                       tipvox, rtipvox,
                       minis, maxes, maxes - minis,
                       [nutvol, nutarea, nutvoxarea],
                       [hull.area, hull.volume, 0]))
    filename = wdst + fname + '_rotation.csv'
    print(filename)
    np.savetxt(filename, align, delimiter=',')
