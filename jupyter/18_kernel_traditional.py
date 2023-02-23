import tifffile as tf
from scipy import ndimage, signal

import numpy as np
import pandas as pd

import os

import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt
import argparse

fs = 15
resol = 0.0759

mborder = ndimage.generate_binary_structure(3, 1).astype(int)
mborder[mborder == 1] = -1
mborder[1,1,1] = -np.sum(mborder) - 1

mborder2 = ndimage.generate_binary_structure(2, 1).astype(int)
mborder2[mborder2 == 1] = -1
mborder2[1,1] = -np.sum(mborder2) - 1

parser = argparse.ArgumentParser(description='Normalize density values of walnuts')

parser.add_argument('src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')
parser.add_argument('dst', metavar='clean_img_dst', type=str, help='path to store clean images')
parser.add_argument('bname', metavar='scan_id', type=str, help='walnut batch scan id')
args = parser.parse_args()

src = args.src
dst = args.dst
bname = args.bname

# src = '../hpcc/'
# dst = '../hpcc/kernel/'

wsrc = src + 'clean/'
tsrc = src + 'watershed/'
rsrc = src + 'rotated/'

walnut_files = sorted(glob.glob(wsrc + bname + '/*.tif'))

wdst = dst + bname + '/'
if not os.path.isdir(wdst):
    os.makedirs(wdst)

for widx in range(len(walnut_files)):
    pa, fi = os.path.split(walnut_files[widx])
    fname = os.path.splitext(fi)[0]
    print(fname)

    tissuefiles = tsrc + bname + '/' + fname + '_'
    meat = tf.imread(tissuefiles + 'meat.tif')

    bimg = meat.copy().astype(int)
    bimg[bimg > 0]  = 1
    border = ndimage.convolve(bimg, mborder, mode='constant', cval=0)
    border[border < 0] = 0
    border[border > 0] = 1

    filename = rsrc + bname + '/' + fname + '_rotation.csv'
    data = np.loadtxt(filename, delimiter=',')

    wmean = data[0]
    rotxyz = data[10:13]
    rotX, _, _ = data[13]
    tipvox = data[14]
    rtipvox = data[15]

    tipvox[0] *= 1 + (-2*rotX)

    datapoints = np.asarray(np.nonzero(border))
    datapoints = datapoints - wmean.reshape(-1,1)
    datapoints *= resol

    mcoords = rotxyz @ datapoints

    kk = 0
    phalfx = mcoords[:, (mcoords[0] >  kk)]
    nhalfx = mcoords[:, mcoords[0] < -kk]

    kk = 1
    phalfy = mcoords[:, mcoords[1] >  kk]
    nhalfy = mcoords[:, mcoords[1] < -kk]

    kk = 1
    phalfz = mcoords[:, (mcoords[2] >  kk)]
    nhalfz = mcoords[:, mcoords[2] < -kk]

    kk = 1
    halfc = mcoords[:, np.abs(mcoords[1]) <  kk]


    H,x,y = np.histogram2d(phalfx[1], phalfx[2], 100, density=False)
    H[H > 0] = 1

    bimg = H.copy().astype(int)
    bimg[bimg > 0]  = 1
    border2 = ndimage.convolve(bimg, mborder2, mode='constant', cval=0)
    border2[border2 < 0] = 0
    border2[border2 > 0] = 1

    c = np.asarray(np.nonzero(border2))
    contour = np.vstack((y[c[1]], x[c[0]]))
    arg = np.argsort(np.sqrt(np.sum(contour**2, axis=0)))[0]

    thetax = np.angle(contour[0,arg] + contour[1,arg]*1j) - np.pi*.5
    rotx = np.array([[np.cos(thetax), -np.sin(thetax)],
                     [np.sin(thetax),  np.cos(thetax)]])
    rcontour = rotx.T @ contour

    rotx3d = np.array([[1,0,0],
                     [0, np.cos(thetax), -np.sin(thetax)],
                     [0, np.sin(thetax),  np.cos(thetax)]])

    rcoords = rotx3d@mcoords

    kk = 0
    rhalfx = rcoords[:, rcoords[0] >  kk]

    kk = 1
    rphalfy = rcoords[:, rcoords[1] >  kk]
    rnhalfy = rcoords[:, rcoords[1] < -kk]

    kk = 1
    rphalfz = rcoords[:, rcoords[2] >  kk]
    rnhalfz = rcoords[:, rcoords[2] < -kk]

    kk = 1
    rhalfc = rcoords[:, np.abs(rcoords[1]) <  kk]

    if np.abs(np.sin(thetax)) > 0.05:
        filename = wdst + fname + '_kernel_rotation.csv'
        print(filename)
        np.savetxt(filename, rotx3d, delimiter=',')

    ################################

    kminis = np.min(mcoords, axis=1)
    kmaxes = np.max(mcoords, axis=1)

    c,b,a = np.sort(kmaxes - kminis)*.5

    krumbein = np.cbrt(b*c/(a*a))
    corey = c/np.sqrt(a*b)
    sneed = np.cbrt(c*c/(a*b))
    janke = c/np.sqrt((a**2 + b**2 + c**2)/3)

    kk = 1
    aux = .75
    mask = (np.abs(mcoords[1]) <= kk) & (np.abs(mcoords[2]) <= kk)
    spine = mcoords[:, mask]
    outer = mcoords[:, ~mask]

    mask = (mcoords[0] < -kk) & (mcoords[1] > kk) & (mcoords[2] > kk)
    foo = mcoords.copy() ; foo[:, ~mask] = 0
    bar = np.abs(foo[0]) - aux*(np.abs(foo[1]))
    leg0x, leg0y, leg0z = mcoords[:, np.argmax(bar)]

    mask = (mcoords[0] < -kk) & (mcoords[1] < -kk) & (mcoords[2] > kk)
    foo = mcoords.copy() ; foo[:, ~mask] = 0
    bar = np.abs(foo[0]) - aux*(np.abs(foo[1]))
    leg1x, leg1y, leg1z = mcoords[:, np.argmax(bar)]

    mask = (mcoords[0] < -kk) & (mcoords[1] < -kk) & (mcoords[2] < -kk)
    foo = mcoords.copy() ; foo[:, ~mask] = 0
    bar = np.abs(foo[0]) - aux*(np.abs(foo[1]))
    leg2x, leg2y, leg2z = mcoords[:, np.argmax(bar)]

    mask = (mcoords[0] < -kk) & (mcoords[1] > kk) & (mcoords[2] < -kk)
    foo = mcoords.copy() ; foo[:, ~mask] = 0
    bar = np.abs(foo[0]) - aux*(np.abs(foo[1]))
    leg3x, leg3y, leg3z = mcoords[:, np.argmax(bar)]

    legx = np.asarray([leg0x, leg1x, leg2x, leg3x])
    legy = np.asarray([leg0y, leg1y, leg2y, leg3y])
    legz = np.asarray([leg0z, leg1z, leg2z, leg3z])

    arm0y = np.max(mcoords[0,(mcoords[1] > kk) & (mcoords[2] > kk)])
    arm0x, arm0y, arm0z = mcoords[:,mcoords[0] == arm0y].squeeze()

    arm1y = np.max(mcoords[0,(mcoords[1] < -kk) & (mcoords[2] > kk)])
    arm1x, arm1y, arm1z = mcoords[:,mcoords[0] == arm1y].squeeze()

    arm2y = np.max(mcoords[0,(mcoords[1] < -kk) & (mcoords[2] < -kk)])
    arm2x, arm2y, arm2z = mcoords[:,mcoords[0] == arm2y].squeeze()

    arm3y = np.max(mcoords[0,(mcoords[1] > kk) & (mcoords[2] < -kk)])
    arm3x, arm3y, arm3z = mcoords[:,mcoords[0] == arm3y].squeeze()

    armx = np.asarray([arm0x, arm1x, arm2x, arm3x])
    army = np.asarray([arm0y, arm1y, arm2y, arm3y])
    armz = np.asarray([arm0z, arm1z, arm2z, arm3z])

    head = np.max(mcoords[0,(np.abs(mcoords[1]) < kk) & (np.abs(mcoords[2]) < kk)])
    head = mcoords[:,mcoords[0] == head].squeeze()

    kk = 1
    hist, bins = np.histogram(spine[0], bins=np.linspace(kminis[0], kmaxes[0], 1000), density = False)
    hist0 = hist.copy()
    fhist = ndimage.median_filter(hist, size=5, mode='constant', cval=0)
    mhistc = ndimage.minimum_filter1d(fhist, size=3, mode='constant', cval=0)

    peaks, _ = signal.find_peaks(hist0, height=np.max(hist0)*.33, wlen=50, distance=125)
    zeros = bins[peaks+1]
    htop = np.max(zeros)
    foo = np.argmin(np.abs(zeros))
    hbot = zeros[foo]
    hlength = htop - hbot

    atop = hbot
    abot = np.median(legx)
    alength = atop - abot


    # SUMARY OF EVERYTHING


    fs = 30

    fig, ax = plt.subplots(3,4,figsize=(25,15), sharex=False, sharey=False, facecolor='snow')
    ax = np.atleast_1d(ax).ravel(); i = 0

    ax[i].imshow(H, origin='lower', cmap='inferno', extent=[y[0],y[-1],x[-1],x[0]])
    ax[i].set_facecolor('k')
    ax[i].set_title(fname, fontsize=fs)
    i += 1
    ax[i].scatter(phalfx[2], phalfx[1], s=.1, color='y', alpha=.2, marker='.')
    ax[i].scatter(contour[0], contour[1], s=25, zorder=3, c='b')
    ax[i].scatter(contour[0,arg], contour[1,arg], marker='^', s=150, c='r', zorder=10)
    ax[i].set_title('{} pts'.format(contour.shape[1]), fontsize=fs)
    i += 1
    ax[i].scatter(rcontour[0], rcontour[1])
    ax[i].set_title('{:.2f} rad'.format(thetax), fontsize=fs)
    i += 1
    ax[i].scatter(rhalfx[2], rhalfx[1], s=.1, color='y', alpha=.2, marker='.')
    ax[i].scatter(rcontour[0], rcontour[1], s=25, zorder=3, c='b')
    ax[i].set_title('{:.2f} deg'.format(np.rad2deg(thetax)), fontsize=fs)
    i += 1

    ######

    ax[i].scatter(phalfz[1], phalfz[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].scatter(legy[:2], legx[:2], s=150, color='m', marker='*');
    ax[i].scatter(army[:2], armx[:2], s=150, color='r', marker='d');
    ax[i].scatter(head[1], head[0], s=150, color='k', marker='^');

    ax[i].plot(legy[:2], np.repeat(abot,2), lw=7, c='m', alpha=0.5);
    ax[i].set_title('{:.2f} mm'.format(legy[0]-legy[1]), fontsize=fs)
    i += 1
    ax[i].scatter(nhalfz[1], nhalfz[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].scatter(legy[2:], legx[2:], s=150, color='m', marker='*');
    ax[i].scatter(army[2:], armx[2:], s=150, color='r', marker='d');
    ax[i].scatter(head[1], head[0], s=150, color='k', marker='^');
    ax[i].plot(legy[2:], np.repeat(abot,2), lw=7, c='m', alpha=0.5);
    ax[i].set_title('{:.2f} mm'.format(legy[3]-legy[2]), fontsize=fs)
    i += 1
    ax[i].scatter(rphalfz[1], rphalfz[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].set_title('Equancy: {:.2f}'.format(c/a), fontsize=fs)
    i+=1
    ax[i].scatter(rnhalfz[1], rnhalfz[0], s=.1, color='y', alpha=.2, marker='.')
    i += 1

    ######

    ax[i].scatter(phalfy[2], phalfy[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].plot([5,5], [htop, hbot], lw=7, c='r', alpha=.5)
    ax[i].plot([-5,-5], [atop, abot], lw=7, c='b', alpha=0.5)
    ax[i].set_title('{:.2f} mm'.format(hlength), fontsize=fs)
    i += 1
    ax[i].scatter(halfc[2], halfc[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].plot([2,2], [htop, hbot], lw=7, c='r', alpha=.5)
    ax[i].plot([-2,-2], [atop, abot], lw=7, c='b', alpha=0.5)
    ax[i].set_title('{:.2f} mm'.format(alength), fontsize=fs)
    i += 1
    ax[i].scatter(rphalfy[2], rphalfy[0], s=.1, color='y', alpha=.2, marker='.')
    i+=1
    ax[i].scatter(rhalfc[2], rhalfc[0], s=.1, color='y', alpha=.2, marker='.')
    i += 1

    for i in range(len(ax)):
        ax[i].set_aspect('equal', 'datalim');
        ax[i].axvline(0, c='g',alpha=0.75)
        ax[i].axhline(0,c='g',alpha=0.75);
    fig.tight_layout();

    filename = wdst + fname + '_kernel.jpg'
    print(filename)
    plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

