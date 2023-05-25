import tifffile as tf
from scipy import ndimage
import numpy as np

import os
import glob
import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import walnut_utils as wnut

src = '../hpcc/clean/'
dst = '../hpcc/watershed/'

struc1 = ndimage.generate_binary_structure(3, 1)
struc2 = ndimage.generate_binary_structure(3, 2)
struc2 = ndimage.generate_binary_structure(3, 3)

ssrange = [100,150,200,250,300]
tissue = ['meat','shell','vein','air']

Sx = [ np.s_[k,:,:] for k in ssrange ]
Sy = [ np.s_[:,k,:] for k in ssrange ]
Sz = [ np.s_[:,:,k] for k in ssrange ]

fs = 15

walnut_batch = sorted(glob.glob( src + '*/' ))

for bidx in range(17,len(walnut_batch)):
#for bidx in [16]:
    walnut_files = sorted(glob.glob(walnut_batch[bidx] + '*.tif'))
    bname = walnut_batch[bidx].split('/')[-2]

    wdst = dst + bname + '/'
    if not os.path.isdir(wdst):
        os.makedirs(wdst)
    wdstd = wdst + 'diagnostic/'
    if not os.path.isdir(wdstd):
        os.makedirs(wdstd)
    wdsts = wdst + 'snaps/'
    if not os.path.isdir(wdsts):
        os.makedirs(wdsts)
    print(wdsts)

    for widx in range(len(walnut_files)):
        print(walnut_files[widx])
        img = tf.imread(walnut_files[widx])

        pa, fi = os.path.split(walnut_files[widx])
        fname = os.path.splitext(fi)[0]
        print(fname)

        mxval = 2**(8*img.itemsize)-1

        #-------------------------------------------------
        # Meat
        #-------------------------------------------------

        mimg = img.copy()
        eimg = ndimage.binary_erosion(mimg, struc2, 15, border_value=1)
        mimg = eimg*mimg
        mimg = ndimage.uniform_filter(mimg)
        mimg[mimg < 120] = 0

        minim = ndimage.minimum_filter(mimg, size=7)

        minimv = minim.flatten()
        minimv = minimv[minimv > 0]

        meanv = np.mean(minimv)
        sdv = np.std(minimv)

        print('Minimum ', meanv, sdv, meanv-2.25*sdv, sep='\t')

        tminim = minim.copy()
        tminim[tminim < meanv-2.25*sdv] = 0
        tminim[tminim > 230] = 0

        eimg = ndimage.binary_erosion(tminim, struc2, 2, border_value=1)
        eimg = wnut.get_largest_element(eimg)

        dimg = ndimage.binary_dilation(eimg, struc2, 6, border_value=0)
        meat = (dimg*img).copy()
        meat[meat < 100] = 0

        #-------------------------------------------------
        # Shell
        #-------------------------------------------------

        dist = ndimage.distance_transform_cdt(img, 'taxicab')
        bshell = (img > 120) & (dist < 7)

        timg = img.copy()
        timg[timg == 0] = mxval
        timg[meat > 0] = 0

        timg = ndimage.binary_erosion(timg, struc2, 5, border_value=1)
        timg[img < 100] = 0
        eimg = wnut.get_largest_element(timg)

        cshell = img*eimg
        cshell[cshell < 180] = 0
        cshell = wnut.get_largest_element(cshell)
        cshell[cshell > 0] = 1

        dshell = ndimage.binary_dilation(cshell, struc2, 2, border_value=0)
        dshell = ndimage.binary_erosion(dshell, struc2, 2, border_value=1)

        mshell = cshell | bshell
        shell = (img*mshell).copy()

        #-------------------------------------------------
        # Veins
        #-------------------------------------------------

        vimg = img.copy()
        vimg[meat > 0] = 0
        vimg[shell > 0] = 0
        vimg = ndimage.binary_erosion(vimg, struc2, 2, border_value=1)
        vein = img*vimg
        vein[vein < 75] = 0
        vein[vein > 190] = 0

        #-------------------------------------------------
        # Air
        #-------------------------------------------------

        aimg = img.copy()

        aimg[img > 60] = 0
        aimg[meat > 0] = 0
        aimg[shell > 0] = 0
        aimg[vein > 0] = 0
        aimg[aimg > 0] = 1
        air = img*aimg

        #-------------------------------------------------
        # Watershed
        #-------------------------------------------------

        markers = np.copy(img).astype(int)
        markers[markers == 0] = -1
        markers[markers > 0] = 0

        for i,maskt in enumerate([meat,shell,vein,air]):
            markers[maskt > 0] = i+1

        ift = ndimage.watershed_ift(img, markers)

        tissues = [None for x in range(4)]
        for i in range(len(tissues)):
            mask = ift == i+1
            box = img.copy()
            box[~mask] = 0

            tissues[i] = box

        #-------------------------------------------------
        # Save and plot
        #-------------------------------------------------


        for i in range(len(tissue)):
            filename = wdst + fname + '_' + tissue[i] + '.tif'
            tf.imwrite(filename, tissues[i], photometric='minisblack', compress=3)

        for i in range(len(tissue)):
            snaps = wnut.collapse_dimensions(tissues[i])
            wnut.plot_collapse_dimensions(snaps, fname, tissue[i], dst=wdsts, writefig=True, display=False)

        # ****************************************************

        fig, ax = plt.subplots(5,5,figsize=(20,21), sharex=True, sharey=True)
        for j in range(5):
            ss = Sx[j]
            i = 0
            ax[j,i].imshow(img[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
            ax[j,i].set_ylabel('slice {}'.format(ssrange[j]), fontsize=fs)
            i = 1
            ax[j,i].imshow(ift[ss], cmap='magma', origin='lower', vmin=0)
            for i in range(2,5):
                ax[j,i].imshow(tissues[i-2][ss], cmap='inferno', origin='lower', vmax=175, vmin=0)
        fig.tight_layout()
        fig.suptitle(fname, fontsize=fs, color='white');

        filename = wdstd + fname + '_x.jpg'
        plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()

        # ****************************************************

        fig, ax = plt.subplots(5,5,figsize=(20,21), sharex=True, sharey=True)
        for j in range(5):
            ss = Sy[j]
            i = 0
            ax[j,i].imshow(img[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
            ax[j,i].set_ylabel('slice {}'.format(ssrange[j]), fontsize=fs)
            i = 1
            ax[j,i].imshow(ift[ss], cmap='magma', origin='lower', vmin=0)
            for i in range(2,5):
                ax[j,i].imshow(tissues[i-2][ss], cmap='inferno', origin='lower', vmax=175, vmin=0)
        fig.tight_layout()
        fig.suptitle(fname, fontsize=fs, color='white');

        filename = wdstd + fname + '_y.jpg'
        plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()


        # ****************************************************

        fig, ax = plt.subplots(5,5,figsize=(20,21), sharex=True, sharey=True)
        for j in range(5):
            ss = Sz[j]
            i = 0
            ax[j,i].imshow(img[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
            ax[j,i].set_ylabel('slice {}'.format(ssrange[j]), fontsize=fs)
            i = 1
            ax[j,i].imshow(ift[ss], cmap='magma', origin='lower', vmin=0)
            for i in range(2,5):
                ax[j,i].imshow(tissues[i-2][ss], cmap='inferno', origin='lower', vmax=175, vmin=0)
        fig.tight_layout()
        fig.suptitle(fname, fontsize=fs, color='white');

        filename = wdstd + fname + '_z.jpg'
        plt.savefig(filename, dpi=96, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
        plt.close()
