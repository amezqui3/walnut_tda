import tifffile as tf
from scipy import ndimage, signal
import numpy as np
import numpy.polynomial.polynomial as P
import os
import argparse
import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import walnut_utils as wnut

dst = '../clean/'
filesrc = '../clean/2014SBa_R1_T25/diagnostic/'
anchors0 = np.loadtxt(filesrc+'anchors0.csv', dtype='float', delimiter=',')
anchors1 = np.loadtxt(filesrc+'anchors1.csv', dtype='float', delimiter=',')
loc_maxa = np.loadtxt(filesrc+'local_maxima.csv', dtype='int', delimiter=',')
y = np.median(loc_maxa, axis=0)
anchors = np.median(anchors1, axis = 0)

walnut_batch = sorted(glob.glob('../raw/*/'))

#for bidx in range(len(walnut_batch)):
for bidx in [4]:
    walnut_files = sorted(glob.glob(walnut_batch[bidx] + '*.tif'))
    bname = walnut_batch[bidx].split('/')[-2]
    wdst = dst + bname + '/'
    if not os.path.isdir(wdst + 'normalization/'):
        os.makedirs(wdst + 'normalization/')

    for widx in range(len(walnut_files)):
        img = tf.imread(walnut_files[widx])//256
        img = img.astype(np.uint8)

        pa, fi = os.path.split(walnut_files[widx])
        fname = os.path.splitext(fi)[0]
        print(walnut_files[widx])

        hist0,bins = np.histogram(img, bins=2**(img.dtype.itemsize*8),range=(0,2**(img.dtype.itemsize*8)))
        fhist = ndimage.median_filter(hist0, size=5, mode='constant', cval=0)
        cumul = np.cumsum(fhist)

        lw=3
        fig = plt.figure(figsize=(15,5))
        plt.plot(np.log(hist0 +1), lw=lw);
        plt.plot(np.log(fhist+1), lw=lw);

        tot = cumul[-1]
        x = np.zeros(len(anchors))
        for i,a in enumerate(anchors):
            for j,s in enumerate(cumul):
                if s>a*tot:
                    break
            if j > 0:
                x[i] = j-1+(a*tot-cumul[j-1])/(cumul[j]-cumul[j-1])
            else:
                x[i] = 0
        npz = P.Polynomial.fit(x,y,2)

        ################################################

        fig, ax = plt.subplots(figsize=(14,5))
        lw = 3

        ax.axhline(0, c='gray', lw=1)
        ax.axvline(y[0], ls='-.', lw=lw, c='r')
        ax.axvline(y[1], ls='-.', lw=lw, c='b')
        ax.axvline(y[2], ls='-.', lw=lw, c='g')
        ax.axvline(0, c='gray', lw=1)
        ax.axvline(255, c='gray', lw=1)
        ax.axvline(75, c='gray', lw=1)

        ax.plot(bins[:-1], np.log(hist0+1), lw=lw, label = 'original')
        ax.plot(bins[:-1], np.log(fhist+1), lw=lw, label = 'filtered')
        ax.plot(npz(bins[:-1]), np.log(fhist+1), lw=lw, label = 'adjusted', c='m')

        ax.set_xlabel("Density", fontsize=22)
        ax.set_ylabel("log(freq+1)", fontsize=22)
        ax.set_title(bname + ' vs. Reference', fontsize=26)
        ax.tick_params(labelsize=22)
        ax.legend(fontsize=24)
        fig.tight_layout()

        filename = wdst + 'normalization/density_normalizaiton_' + fname + '.jpg'
        plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()

        ################################################

        aimg = img.copy()
        aimg = wnut.normalize_density(aimg, npz)

        thr = 100
        struc = ndimage.generate_binary_structure(img.ndim, 2)
        cimg = aimg.copy()
        cimg[cimg < thr] = 0
        oimg = wnut.fill_component(cimg)
        oimg = wnut.get_largest_element(oimg)
        cimg = (oimg > 0)*aimg


        ################################################

        k = 250
        ss = np.s_[k,:, :]

        fig, ax = plt.subplots(2,2,figsize=(15,15), sharex=True, sharey=True)
        ax = np.atleast_1d(ax).flatten()

        i = 0
        ax[i].imshow(img[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('img')
        i = 1
        ax[i].imshow(aimg[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('normalized')
        i = 2
        ax[i].imshow(oimg[ss], cmap='inferno', origin='lower', vmax=1, vmin=0)
        ax[i].set_title('components filled')
        i = 3
        ax[i].imshow(cimg[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('air removed')

        fig.tight_layout()
        filename = wdst + 'normalization/air_removal_' + fname + '.jpg'
        plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()

        walnut, cero = wnut.clean_zeroes(cimg)
        np.savetxt(wdst + 'normalization/clean_zeroes.csv', cero, fmt='%d', delimiter = ',')

        filename = wdst + fname + '.tif'
        tf.imwrite(filename, walnut, photometric='minisblack', compress=3)

        snaps = wnut.collapse_dimensions(walnut)
        wnut.plot_collapse_dimensions(snaps, fname, 'whole walnut', dst=wdst+'normalization/', writefig=True)

        snaps = wnut.collapse_dimensions_max(walnut)
        wnut.plot_collapse_dimensions(snaps, fname, 'walnut shell', dst=wdst+'normalization/', writefig=True)
