import tifffile as tf
from scipy import ndimage
from scipy import spatial
from scipy import special
import numpy as np
import pandas as pd

import os
import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import walnut_utils as wnut
import argparse

fs = 15
resol = 0.0759

mborder = ndimage.generate_binary_structure(3, 1).astype(int)
mborder[mborder == 1] = -1
mborder[1,1,1] = -np.sum(mborder) - 1

parser = argparse.ArgumentParser(description='Normalize density values of walnuts')

parser.add_argument('src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')
parser.add_argument('dst', metavar='clean_img_dst', type=str, help='path to store clean images')
parser.add_argument('bname', metavar='scan_id', type=str, help='walnut batch scan id')
args = parser.parse_args()

src = args.src
dst = args.dst
bname = args.bname

# src = '../hpcc/'
# dst = '../hpcc/traditional/'

wsrc = src + 'clean/'
tsrc = src + 'watershed/'
rsrc = src + 'rotated/'

walnut_files = sorted(glob.glob(wsrc + bname + '/*.tif'))

wdst = dst + bname + '/'
if not os.path.isdir(wdst):
    os.makedirs(wdst)

wdstp = wdst + 'protrusion/'
if not os.path.isdir(wdstp):
    os.makedirs(wdstp)

for widx in range(len(walnut_files)):
    img = tf.imread(walnut_files[widx])
    pa, fi = os.path.split(walnut_files[widx])
    fname = os.path.splitext(fi)[0]

    tissuefiles = tsrc + bname + '/' + fname + '_'
    air  = tf.imread(tissuefiles + 'air.tif')
    meat = tf.imread(tissuefiles + 'meat.tif')
    shell= tf.imread(tissuefiles + 'shell.tif')
    vein = tf.imread(tissuefiles + 'vein.tif')
    protrusion = tf.imread(tissuefiles + 'protrusion.tif')

    extshell = np.zeros_like(shell, dtype=np.bool)
    extshell[protrusion == 3] = True

    filename = rsrc + bname + '/' + fname + '_rotation.csv'
    data = np.loadtxt(filename, delimiter=',')

    wmean = data[0]
    rotxyz = data[10:13]
    rotX, _, _ = data[13]
    tipvox = data[14]
    rtipvox = data[15]
    feretd = data[18]
    nutvol, nutarea, nutvoxarea = data[19]
    chnutarea, chnutvol, _ = data[20]

    tipvox[0] *= 1 + (-2*rotX)

    chnutaratio = chnutarea/nutarea
    chnutvratio = chnutvol/nutvol

    # ### Other phenotypes

    tvols = np.zeros(4)
    for i,tissue in enumerate([air, meat, shell, vein]):
        tvols[i] = np.sum(tissue > 0)

    tvols = tvols.astype(float)*(resol**3)

    tvolr = tvols/nutvol

    # ### Sphericity

    nutva3d = (nutarea ** 3)/(36*np.pi*nutvol**2)
    nutferet = np.max(feretd)/np.min(feretd)
    nutsphr = np.cbrt(36 * np.pi * nutvol**2)/nutarea
    shellrug = 1./nutsphr

    c,b,a = np.sort(feretd)*.5

    if a == c:
        area = 4*np.pi*a*a
    else:
        t = np.arccos(c/a)
        s = np.arccos(c/b)
        k = np.sin(s)/np.sin(t)

        ellarea = c*c/(a*a)*special.ellipkinc(t, k*k) + np.sin(t)*np.sin(t)*special.ellipeinc(t, k*k)
        ellarea *= 2*np.pi*a*b/np.sin(t)
        ellarea += 2*np.pi*c*c

    ellvolume = 4*np.pi*a*b*c/3

    wadell = np.cbrt(36*np.pi*ellvolume*ellvolume)/ellarea

    krumbein = np.cbrt(b*c/(a*a))
    corey = c/np.sqrt(a*b)
    sneed = np.cbrt(c*c/(a*b))
    janke = c/np.sqrt((a**2 + b**2 + c**2)/3)

    # ### Kernel lobeyness

    bimg = meat.copy().astype(int)
    bimg[bimg > 0]  = 1
    border = ndimage.convolve(bimg, mborder, mode='constant', cval=0)
    border[border < 0] = 0
    kerarea = np.sum(border) * (resol ** 2)
    border[border > 0] = 1
    kervoxarea = np.sum(border)

    datapoints = np.asarray(np.nonzero(border))
    datapoints = datapoints - wmean.reshape(-1,1)
    datapoints *= resol

    mcoords = rotxyz @ datapoints

    kminis = np.min(mcoords, axis=1)
    kmaxes = np.max(mcoords, axis=1)

    khull = spatial.ConvexHull(mcoords.T)
    kerlob = khull.area/kerarea
    chkervratio = khull.volume/tvols[1]

    # ### Shell thickness

    thickness = wnut.object_thickness(extshell, resol)

    # ## Protruding shell

    shellvols, _ = np.histogram(protrusion, [2, 4, 6, 10], range=(0,10))

    # # Data saving

    tradpheno = np.hstack((feretd,
                           nutvol,
                           nutva3d, # rugosity ** 3
                           nutferet,
                           1./nutferet, # equancy
                           nutarea,
                           nutsphr, # wadell
                           chnutarea,
                           chnutvol,
                           chnutaratio,
                           chnutvratio,
                           1./chnutaratio,
                           1./chnutvratio,
                           krumbein,
                           corey,
                           sneed,
                           janke,
                           wadell,
                           tvols,
                           tvolr,
                           shellrug, # 1/wadell
                           thickness,
                           shellvols/np.sum(shellvols),
                           shellvols*resol**3,
                           kmaxes - kminis,
                           kerarea,
                           khull.volume,
                           khull.area,
                           kerlob,
                           1./kerlob,
                           chkervratio,
                           1./chkervratio))

    filename = wdst + fname + '_trad.csv'
    foo = pd.DataFrame([bname, fname.split('_')[-1], *tradpheno]).T
    foo.to_csv(filename, header=False, index=False)
