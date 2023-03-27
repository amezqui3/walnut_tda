import tifffile as tf
from scipy import ndimage, signal, spatial

import numpy as np
import pandas as pd
import os

import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt
import argparse
from shapely.geometry import LineString, Polygon
from shapely.ops import polylabel

ang = np.linspace(0, 2*np.pi, 200, endpoint=True)
circle = np.vstack((np.sin(ang), np.cos(ang)))

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
    kervol = np.sum(bimg)*(resol**3)
    border = ndimage.convolve(bimg, mborder, mode='constant', cval=0)
    border[border < 0] = 0
    kerarea = np.sum(border) * (resol ** 2)
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
    khull = spatial.ConvexHull(mcoords.T)
    charearatio = kerarea/khull.area
    chivolratio = khull.volume/kervol

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
    halfy = mcoords[:, np.abs(mcoords[1]) <  kk]
    halfz = mcoords[:, np.abs(mcoords[2]) <  kk]

    kminis = np.min(mcoords, axis=1)
    kmaxes = np.max(mcoords, axis=1)

    c,b,a = np.sort(kmaxes - kminis)*.5

    krumbein = np.cbrt(b*c/(a*a))
    corey = c/np.sqrt(a*b)
    sneed = np.cbrt(c*c/(a*b))
    janke = c/np.sqrt((a**2 + b**2 + c**2)/3)

    kk = 1
    mask = (np.abs(mcoords[1]) <= kk) & (np.abs(mcoords[2]) <= kk)
    spine = mcoords[:, mask]
    outer = mcoords[:, ~mask]

    kk = 1
    aux = 1
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

    # central bit

    kk = 1
    hist, bins = np.histogram(spine[0], bins=np.linspace(kminis[0], kmaxes[0], 1000), density = False)
    hist0 = hist.copy()
    fhist = ndimage.median_filter(hist, size=5, mode='constant', cval=0)
    mhistc = ndimage.minimum_filter1d(fhist, size=3, mode='constant', cval=0)

    peaks, _ = signal.find_peaks(mhistc, height=np.max(mhistc)*.33, prominence=np.max(mhistc)*.33, wlen=50)
    if len(peaks) < 2:
        peaks, _ = signal.find_peaks(hist0, height=np.max(hist0)*.33, prominence=np.max(hist0)*.33, wlen=50)
    zeros = bins[peaks+1]
    htop = np.max(zeros)
    foo = np.argmin(np.abs(zeros))
    hbot = zeros[foo]
    if hbot < -4:
        hbot = 0
    if htop < hbot + 7:
        htop = bins[np.nonzero(mhistc)[0][-1]]
    hlength = htop - hbot

    # lower arch height
    atop = hbot
    abot = np.median(legx)
    alength = atop - abot
    print(atop, abot, alength, sep='\t')

    linspace = np.linspace(atop-1, abot+1, 25)
    eps = (linspace[0] - linspace[1])/2
    granularity = 100
    bins = np.linspace(-np.pi, np.pi, granularity+1)

    surface = np.zeros(len(linspace))
    volume = np.zeros_like(surface)
    distantpoles = np.zeros_like(volume)
    distantpolev = np.zeros_like(volume)

    fs = 14
    fig, ax = plt.subplots(5,5,figsize=(16,16), sharex=True, sharey=True, facecolor='snow')
    ax = np.atleast_1d(ax).ravel(); i = 0

    for i in range(len(linspace)):
        kk = linspace[i]
        nhalfx = mcoords[:, (mcoords[0] < kk + eps) & (mcoords[0] > kk - eps)]
        #nhalfx = nhalfx - np.mean(nhalfx, axis=1).reshape(-1,1)

        data = nhalfx[1:]
        angles = np.angle(data[0] + data[1]*1j)

        binning = np.digitize(angles, bins)

        trace = np.zeros((granularity, 2))
        for j in range(granularity):
            if np.sum(binning == j + 1) > 1:
                subdata = data[:, binning == j + 1]
                arg = np.argmin(np.sqrt(np.sum(subdata**2, axis=0)))
                trace[j] = subdata[:, arg]
        trace = trace[np.all(trace, axis=1),:].T
        radii = np.sqrt(np.sum(trace**2, axis=0))
        mr = np.quantile(radii, 0.8)

        trace = trace[:, radii < mr]
        polygon = Polygon(zip(*trace))
        center = polylabel(polygon, tolerance=1e-1)
        rr = np.min(np.sqrt((trace[0] - center.x)**2 + (trace[1] - center.y)**2))

        surface[i] = polygon.length*eps
        volume[i] = polygon.area*eps
        distantpoles[i] = 4*np.pi*rr**2
        distantpolev[i] = 4/3*np.pi*rr**3

        ax[i].set_title('A: {:.1f}, r: {:.1f}'.format(polygon.area, rr), fontsize=fs)
        ax[i].scatter(nhalfx[1], nhalfx[2], s=0.5, color='y', alpha=1, marker='.');
        ax[i].plot(trace[0], trace[1], c='r', lw=2)
        ax[i].plot([trace[0,0],trace[0,-1]], [trace[1,0],trace[1,-1]], c='r', lw=2)
        ax[i].scatter([center.x], [center.y], s=50, c='k')
        ax[i].plot(rr*circle[0] + center.x, rr*circle[1]+center.y, lw=2, c='k', zorder=5)

    for i in range(len(ax)):
        ax[i].set_aspect('equal');
        ax[i].axvline(0, c='g',alpha=0.75)
        ax[i].axhline(0,c='g',alpha=0.75);

    fig.suptitle(fname, fontsize=fs+10)
    fig.tight_layout();
    filename = wdst + fname + '_cavity.jpg'
    print(filename)
    plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

    vol = np.sum(volume)
    suf = np.sum(surface)
    fs = 30

    fig, ax = plt.subplots(2,4,figsize=(20,10), sharex=False, sharey=False, facecolor='snow')
    ax = np.atleast_1d(ax).ravel(); i = 0

    ax[i].scatter(phalfx[1], phalfx[2], s=.1, color='y', alpha=.2, marker='.')
    ax[i].set_title(fname, fontsize=fs)
    i += 1
    ax[i].scatter(phalfz[1], phalfz[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].axhline(np.mean([atop,abot]), c='b', lw=5, ls='--', alpha=.5)
    ax[i].set_title('Equancy: {:.2f}'.format(c/a), fontsize=fs)
    i += 1
    ax[i].scatter(phalfy[2], phalfy[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].axhline(np.mean([atop,abot]), c='b', lw=5, ls='--', alpha=.5)
    ax[i].set_title('Janke: {:.2f}'.format(janke), fontsize=fs)
    i += 1
    ax[i].scatter(halfy[2], halfy[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].plot([2,2], [htop, hbot], lw=7, c='r', alpha=.5)
    ax[i].plot([-2,-2], [atop, abot], lw=7, c='b', alpha=0.5)
    ax[i].set_title('red: {:.2f} mm'.format(hlength), fontsize=fs)
    i += 1

    ######

    ax[i].scatter(nhalfx[1], nhalfx[2], s=.1, color='y', alpha=.2, marker='.')
    ax[i].scatter(trace[0], trace[1], s=10, zorder=4, c='r', marker='*')
    ax[i].set_title('X: {:.2f}'.format(vol/suf), fontsize=fs)
    i += 1
    ax[i].scatter(nhalfz[1], nhalfz[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].axhline(np.mean([atop,abot]), c='b', lw=5, ls='--', alpha=.5)
    ax[i].set_title('S: {:.0f} mm2'.format(suf), fontsize=fs)
    i += 1
    ax[i].scatter(nhalfy[2], nhalfy[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].axhline(np.mean([atop,abot]), c='b', lw=5, ls='--', alpha=.5)
    ax[i].set_title('V: {:.0f} mm3'.format(vol), fontsize=fs)
    i += 1
    ax[i].scatter(halfz[1], halfz[0], s=.1, color='y', alpha=.2, marker='.')
    ax[i].plot([2,2], [htop, hbot], lw=7, c='r', alpha=.5)
    ax[i].plot([-2,-2], [atop, abot], lw=7, c='b', alpha=0.5)
    ax[i].set_title('blue: {:.2f} mm'.format(alength), fontsize=fs)

    for i in range(len(ax)):
        ax[i].set_aspect('equal', 'datalim');
        ax[i].axvline(0, c='g',alpha=0.75)
        ax[i].axhline(0,c='g',alpha=0.75);
    fig.tight_layout();

    filename = wdst + fname + '_kernel.jpg'
    print(filename)
    plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

    kerpheno = [kerarea, kervol, charearatio, chivolratio, krumbein,
            corey, sneed, janke, c/a, suf, vol, hlength, alength]
    filename = wdst + fname + '_kernel.csv'
    foo = pd.DataFrame([bname, fname.split('_')[-1], *kerpheno]).T
    foo.to_csv(filename, header=False, index=False)
