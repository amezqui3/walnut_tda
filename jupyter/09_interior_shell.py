import tifffile as tf
from scipy import ndimage
import numpy as np

import os
import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import walnut_utils as wnut
import argparse

fs = 15
resol = 0.0759
c = 'w'
R = 5
PAD = 30

struc1 = ndimage.generate_binary_structure(3, 1)
struc2 = ndimage.generate_binary_structure(3, 2)
struc3 = ndimage.generate_binary_structure(3, 3)

ssrange = [50,220,240,260,280,-75]
rrrange = np.linspace(0, 2*np.pi, 100)
xcirc, ycirc = R*np.cos(rrrange), R*np.sin(rrrange)

Sx = [ np.s_[k,:,:] for k in ssrange ]
Sy = [ np.s_[:,k,:] for k in ssrange ]
Sz = [ np.s_[:,:,k] for k in ssrange ]
Sxyz = Sx + Sy + Sz

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
# dst = '../hpcc/watershed/'

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

for widx in range(7,len(walnut_files)):
    img = tf.imread(walnut_files[widx])
    pa, fi = os.path.split(walnut_files[widx])
    fname = os.path.splitext(fi)[0]
    print(fname)

    shell= tf.imread(tsrc + bname + '/' + fname + '_shell.tif')

    distL1 = ndimage.distance_transform_cdt(shell, 'taxicab')
    K = np.max(distL1)

    filename = rsrc + bname + '/' + fname + '_rotation.csv'
    data = np.loadtxt(filename, delimiter=',')

    wmean = data[0]
    rotxyz = data[10:13]

    bimg = img.copy().astype(int)
    bimg[bimg > 0]  = 1
    border = ndimage.convolve(bimg, mborder, mode='constant', cval=0)
    border[border < 0] = 0
    nutarea = np.sum(border) * (resol ** 2)
    border[border > 0] = 1
    orig_datapoints = np.asarray(np.nonzero(border))
    mcoords = rotxyz @ ((orig_datapoints - wmean.reshape(-1,1))*resol)

    dshell = np.pad(shell, PAD)
    dshell = ndimage.binary_dilation(dshell, struc2, PAD-5, border_value=0)

    eshell = ndimage.binary_erosion(dshell, struc2, PAD-5, border_value=1)
    eshell = eshell[PAD:-PAD, PAD:-PAD, PAD:-PAD]
    fshell = wnut.fill_component(eshell)
    efill0 = ndimage.binary_erosion(fshell, struc1, K)
    efill = ndimage.binary_erosion(fshell, struc1, 2*K+10)

    extra = efill0*shell

    test = np.zeros_like(img)
    test[shell > 0] += 3
    test[extra > 0] += 2

    ishell = efill*shell
    if np.sum(ishell > 0) == 0:
        efill = ndimage.binary_erosion(fshell, struc1, 2*K)
        ishell = efill*shell
    ishell, labels, where = wnut.get_largest_element(ishell,1e-3, outlabels=True)
    coms = ndimage.center_of_mass(ishell, labels, where)
    coms = np.asarray(coms).T

    mcoms = rotxyz @ ((coms - wmean.reshape(-1,1))*resol)
    dmcoms = np.sqrt((mcoms[1]**2 + mcoms[2]**2))
    mwhere = np.where( dmcoms < R)[0]
    nwhere = np.where( dmcoms >= R)[0]

    if(len(mwhere) > 0):
        fig, ax = wnut.plot_3Dprojections(mcoords, fname);
        ax[0].plot(xcirc, ycirc, c='gray', lw=3, ls='--')
        ax[1].axvline(R, c='gray', lw=3, ls='--'); ax[1].axvline(-R, c='gray', lw=3, ls='--', zorder=1)
        ax[2].axvline(R, c='gray', lw=3, ls='--'); ax[2].axvline(-R, c='gray', lw=3, ls='--', zorder=1)

        for i in range(len(mwhere)):
            ax[0].scatter(mcoms[2,mwhere[i]], mcoms[1,mwhere[i]], marker='*', s=50, c='b', alpha=1, zorder=3)
            ax[1].scatter(mcoms[2,mwhere[i]], mcoms[0,mwhere[i]], marker='*', s=50, c='b', alpha=1, zorder=3)
            ax[2].scatter(mcoms[1,mwhere[i]], mcoms[0,mwhere[i]], marker='*', s=50, c='b', alpha=1, zorder=3);
        for i in range(len(nwhere)):
            ax[0].scatter(mcoms[2,nwhere[i]], mcoms[1,nwhere[i]], marker='^', s=35, c='r', alpha=1, zorder=2)
            ax[1].scatter(mcoms[2,nwhere[i]], mcoms[0,nwhere[i]], marker='^', s=35, c='r', alpha=1, zorder=2)
            ax[2].scatter(mcoms[1,nwhere[i]], mcoms[0,nwhere[i]], marker='^', s=35, c='r', alpha=1, zorder=2);
        filename = wdstp + fname + '_coords'
        plt.savefig(filename + '.jpg', dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

        imask = np.zeros_like(img).astype(np.bool)
        for i in where[mwhere]:
            imask[labels == i] = True

        imask = ndimage.binary_dilation(imask, struc3, int(np.ceil(K*.1)))
        imask = ndimage.binary_dilation(imask, struc1, int(np.ceil(K*.5)))
        protruding = shell*imask
        test[protruding > 0] += 4

        fig, ax = plt.subplots(3,len(ssrange),figsize=(14,8), sharex=False, sharey=False, facecolor='k')
        for i in range(3):
            for j in range(len(ssrange)):
                ss = Sxyz[len(ssrange)*i + j]
                ax[i,j].imshow(imask[ss] + shell[ss], cmap='magma', origin='lower', vmax=2, vmin=0)
                ax[i,j].get_xaxis().set_ticks([])
                ax[i,j].get_yaxis().set_ticks([])
                for spine in ax[i,j].spines.values():
                    spine.set_visible(False)
        for i,x in enumerate(['X','Y','Z']):
            ax[i,0].set_ylabel(x, fontsize=fs, color=c, rotation='horizontal', ha='center', va='center')
        for j in range(len(ssrange)-1):
            ax[-1,j].set_xlabel(ssrange[j], fontsize=fs, color=c)
        xlabel = '/'.join((np.asarray(img.shape) + ssrange[-1]).astype(str))
        ax[-1,-1].set_xlabel(xlabel, color=c, fontsize=fs)

        fig.suptitle(fname, fontsize=fs+5, color=c);
        fig.tight_layout()

        filename = wdstp + fname + '_detection'
        plt.savefig(filename + '.jpg', dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

    fig, ax = plt.subplots(3,len(ssrange),figsize=(14,8), sharex=False, sharey=False, facecolor='k')
    for i in range(3):
        for j in range(len(ssrange)):
            ss = Sxyz[len(ssrange)*i + j]
            ax[i,j].imshow(test[ss], cmap='magma', origin='lower', vmax=7, vmin=0)
            ax[i,j].get_xaxis().set_ticks([])
            ax[i,j].get_yaxis().set_ticks([])
            for spine in ax[i,j].spines.values():
                spine.set_visible(False)
    for i,x in enumerate(['X','Y','Z']):
        ax[i,0].set_ylabel(x, fontsize=fs, color=c, rotation='horizontal', ha='center', va='center')
    for j in range(len(ssrange)-1):
        ax[-1,j].set_xlabel(ssrange[j], fontsize=fs, color=c)
    xlabel = '/'.join((np.asarray(img.shape) + ssrange[-1]).astype(str))
    ax[-1,-1].set_xlabel(xlabel, color=c, fontsize=fs)

    fig.suptitle(fname, fontsize=fs+5, color=c);
    fig.tight_layout()

    filename = wdstp + fname + '_protrusion'
    plt.savefig(filename + '.jpg', dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})

    snaps = wnut.collapse_dimensions_max(test)
    fig, ax = wnut.plot_collapse_dimensions(snaps, fname, 'cavity', dst=wdstp, writefig=False, display=True)
    for i in range(len(ax)):
        ax[i].axis('off')
    filename = wdstp + fname + '_cavity.jpg'
    plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight');

    filename = wdst + fname + '_protrusion.tif'
    tf.imwrite(filename, test, photometric='minisblack', compress=3)

