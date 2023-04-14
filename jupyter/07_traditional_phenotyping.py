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

parser.add_argument('raw_src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')
parser.add_argument('src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')
parser.add_argument('dst', metavar='clean_img_dst', type=str, help='path to store clean images')
parser.add_argument('bname', metavar='scan_id', type=str, help='walnut batch scan id')
args = parser.parse_args()

raw_src = args.raw_src
src = args.src
dst = args.dst
bname = args.bname

# raw_src = '../raw/'
# src = '../hpcc/'
# dst = '../hpcc/traditional/'

wsrc = src + 'clean/'
tsrc = src + 'watershed/'
rsrc = src + 'rotated/'
ksrc = src + 'kernel/'

walnut_files = sorted(glob.glob(wsrc + bname + '/*.tif'))

wdst = dst + bname + '/'
if not os.path.isdir(wdst):
    os.makedirs(wdst)

#for widx in range(len(walnut_files)):
for widx in range(2,len(walnut_files)):

    pa, fi = os.path.split(walnut_files[widx])
    fname = os.path.splitext(fi)[0]
    lname = fname[-3:]

    print(fname)

    img = tf.imread(wsrc + bname + '/' + fname + '.tif')
    tissuefiles = tsrc + bname + '/' + fname + '_'
    air  = tf.imread(tissuefiles + 'air.tif')
    meat = tf.imread(tissuefiles + 'meat.tif')
    shell= tf.imread(tissuefiles + 'shell.tif')
    vein = tf.imread(tissuefiles + 'vein.tif')
    protrusion = tf.imread(tissuefiles + 'protrusion.tif')
    extshell = np.zeros_like(shell, dtype=np.bool)
    extshell[protrusion == 3] = True

    ## Relative densities

    raw_file = raw_src + bname + '/' + lname + '.tif'
    print(raw_file)
    raw = tf.imread(raw_file)

    zerofile = wsrc + bname + '/normalization/clean_zeroes' + lname + '.csv'
    cero = np.loadtxt(zerofile, dtype=int)
    raw = raw[cero[1]:cero[4], cero[2]:cero[5], cero[0]:cero[3]]

    if not raw.shape == meat.shape:
        raw = raw.swapaxes(0,1)

    if not raw.shape == meat.shape:
        raw = raw.swapaxes(0,1)
        raw = raw.swapaxes(0,2)

    if not raw.shape == meat.shape:
        raw = raw.swapaxes(0,2)
        raw = raw.swapaxes(1,2)

    raw_meat = raw.copy()
    raw_meat[meat == 0] = 0

    raw_shell = raw.copy()
    raw_shell[shell == 0] = 0

    raw_vein = raw.copy()
    raw_vein[vein == 0] = 0

    rho_meat = np.mean(raw_meat[raw_meat > 0])
    rho_shell = np.mean(raw_shell[raw_shell > 0])
    rho_vein = np.mean(raw_vein[raw_vein > 0])

    rho_mvs = rho_meat/rho_shell
    rho_vvs = rho_vein/rho_shell
    rho_vvm = rho_vein/rho_meat

    ## Load rotations and phenotypes

    filename = rsrc + bname + '/' + fname + '_rotation.csv'
    data = np.loadtxt(filename, delimiter=',')

    wmean = data[0]
    rotxyz = data[10:13]
    rotX, _, _ = data[13]
    feretd = data[18]
    nutvol, nutarea, nutvoxarea = data[19]
    chnutarea, chnutvol, _ = data[20]

    chnutaratio = chnutarea/nutarea
    chnutvratio = chnutvol/nutvol

    ## Other phenotypes

    tvols = np.zeros(4)
    for i,tissue in enumerate([air, meat, shell, vein]):
        tvols[i] = np.sum(tissue > 0)

    tvols = tvols.astype(float)*(resol**3)
    tvolr = tvols/nutvol

    ## Sphericity

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
    wwadell = np.cbrt(36*np.pi*ellvolume*ellvolume)/ellarea

    wkrumbein = np.cbrt(b*c/(a*a))
    wcorey = c/np.sqrt(a*b)
    wsneed = np.cbrt(c*c/(a*b))
    wjanke = c/np.sqrt((a**2 + b**2 + c**2)/3)
    wequancy = c/a

    ## Kernel phenotypes

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

    filename = ksrc + bname + '/' + fname + '_kernel.csv'
    kpheno = pd.read_csv(filename, header=None, dtype={1:str}).values[0][2:].astype(float)

    # 0   kerarea
    # 1   kervol
    # 2   charearatio
    # 3   chivolratio
    # 4   krumbein
    # 5   corey
    # 6   sneed
    # 7   janke
    # 8   c/a
    # 9   suf
    # 10  vol
    # 11  hlength
    # 12  alength

    kerSR = kpheno[9]/kpheno[0]
    kerVR = kpheno[10]/khull.volume
    kersphr = np.cbrt(36 * np.pi * kpheno[1]**2)/kpheno[0]

    ## Shell thickness

    thickness = wnut.object_thickness(extshell, resol, NNN=4, K=5)[0]

    shellvols, _ = np.histogram(protrusion, [2, 4, 6, 10], range=(0,10))
    shellvols = np.asarray([shellvols[0], shellvols[1]+shellvols[2]])

    shellvols/np.sum(shellvols)
    shellvols*resol**3

    ## Data saving

    tradpheno = np.hstack((feretd,
                           nutvol,
                           nutva3d, # rugosity ** 3
                           nutferet,
                           nutarea,
                           nutsphr, # wadell
                           chnutarea,
                           chnutvol,
                           chnutaratio,
                           1./chnutvratio,
                           wkrumbein,
                           wsneed,
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
                           kpheno[2:4],
                           kpheno[9:],
                           kersphr,
                           kerSR,
                           kerVR,
                           rho_mvs, rho_vvs, rho_vvm))

    filename = wdst + fname + '_trad.csv'
    foo = pd.DataFrame([bname, fname.split('_')[-1], *tradpheno]).T

    print(filename)

    foo.to_csv(filename, header=False, index=False)

