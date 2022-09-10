import tifffile as tf
import scipy.ndimage as ndimage
import numpy as np
import os
import glob
from matplotlib import pyplot as plt
import unionfind as UF
import walnut_utils as wnut

dst = '../clean/'
xaxis = np.arange(3)

walnut_batch = sorted(glob.glob('../raw/*/'))

for bidx in range(len(walnut_batch)):
#for bidx in range(1):
    walnut_files = sorted(glob.glob(walnut_batch[bidx] + '*.tif'))
    bname = walnut_batch[bidx].split('/')[-2]

    wdst = dst + bname + '/'
    if not os.path.isdir(wdst):
        os.makedirs(wdst)
    wdst += 'diagnostic/'
    if not os.path.isdir(wdst):
        os.makedirs(wdst)

    LocMax = []
    Anchor0 = []
    Anchor1 = []

    for widx in range(len(walnut_files)):
    #for widx in range(1):

        pa, fi = os.path.split(walnut_files[widx])
        fname = os.path.splitext(fi)[0]

        k = 250
        ss = np.s_[:, k, :]

        print(walnut_files[widx])
        img = tf.imread(walnut_files[widx])//256
        img = img.astype(np.uint8)

        hist0, bins = np.histogram(img, bins=2**(img.dtype.itemsize*8),range=(0,2**(img.dtype.itemsize*8)))
        fhist = ndimage.median_filter(hist0, size=5, mode='constant', cval=0)

        cumul0 = np.cumsum(hist0)
        cumul1 = np.cumsum(fhist)

        tot0 = cumul0[-1]
        tot1 = cumul1[-1]

        pers = sorted(UF.persistence(fhist),reverse=True)

        p0,p1,p2 = pers[:3]
        x0,x1,x2 = p0[2],p1[2],p2[2]
        anchors0 = [cumul0[x0]/tot0,cumul0[x1]/tot0,cumul0[x2]/tot0]
        anchors1 = [cumul1[x0]/tot1,cumul1[x1]/tot1,cumul1[x2]/tot1]

        LocMax.append([x0,x1,x2])
        Anchor0.append(anchors0)
        Anchor1.append(anchors1)

        ################################################

        fig, ax = plt.subplots(1,1,figsize=(20,9), sharex=True)
        ax = np.atleast_1d(ax).flatten()
        lw = 4

        i = 0
        ax[i].axvline(x0, ls='-.', lw=lw, c='r', label=int(x0));
        ax[i].axvline(x1, ls='-.', lw=lw, c='b', label=int(x1));
        ax[i].axvline(x2, ls='-.', lw=lw, c='g', label=int(x2))
        ax[i].set_title('img : pers')
        ax[i].plot(bins[:-1],np.log(hist0+1), lw=lw-1)
        ax[i].plot(bins[:-1],np.log(fhist+1), lw=lw-1)


        for i in range(len(ax)):
            ax[i].set_xlabel("Density", fontsize=22)
            ax[i].set_ylabel("log(freq+1)", fontsize=22)
            ax[i].tick_params(labelsize=22)
            ax[i].legend(fontsize=12)
        fig.tight_layout()

        filename = wdst + 'density_distribution_' + fname + '.jpg'
        plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()

        ################################################

        cimg0 = img[ss].copy()
        cimg0[cimg0 < x0] = 0

        cimg1 = img[ss].copy()
        cimg1[cimg1 < x1] = 0

        cimg2 = img[ss].copy()
        cimg2[cimg2 < x2] = 0

        fig, ax = plt.subplots(2,2,figsize=(14,18), sharex=True, sharey=True)
        ax = np.atleast_1d(ax).flatten()

        i = 0
        ax[i].imshow(img[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('img')
        i = 1
        ax[i].imshow(cimg0, cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('img[img < {}] = 0'.format(int(x0)))
        i = 2
        ax[i].imshow(cimg1, cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('img[img < {}] = 0'.format(int(x1)))
        i = 3
        ax[i].imshow(cimg2, cmap='inferno', origin='lower', vmax=255, vmin=0)
        ax[i].set_title('img[img < {}] = 0'.format(int(x2)))

        fig.tight_layout()

        filename = wdst + 'walnuts_' + fname + '.jpg'
        plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()

    LocMax = np.asarray(LocMax)
    Anchor0 = np.asarray(Anchor0)
    Anchor1 = np.asarray(Anchor1)


    np.savetxt(wdst + 'local_maxima.csv', LocMax, fmt='%d', delimiter = ',')
    np.savetxt(wdst + 'anchors0.csv', Anchor0, fmt='%.18e', delimiter = ',')
    np.savetxt(wdst + 'anchors1.csv', Anchor1, fmt='%.18e', delimiter = ',')

    fig, ax = plt.subplots(3,1,figsize=(14,10), sharex=True)
    ax = np.atleast_1d(ax).flatten()

    for j in range(len(LocMax)):
        ax[0].scatter(xaxis, LocMax[j], marker='${}$'.format(j+1), s=100, c=['r','g','b'])
        ax[1].scatter(xaxis, Anchor0[j], marker='${}$'.format(j+1), s=100, c=['r','g','b'])
        ax[2].scatter(xaxis, Anchor1[j], marker='${}$'.format(j+1), s=100, c=['r','g','b'])

    ax[0].set_title('Loc Max', fontsize=20)
    ax[1].set_title('Anchor0 points', fontsize=20)
    ax[2].set_title('Anchor1 points', fontsize=20)
    fig.tight_layout()

    filename = wdst + 'anchor_points_' + bname + '.jpg'
    plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
    plt.close()
