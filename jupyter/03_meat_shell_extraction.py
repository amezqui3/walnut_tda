import tifffile as tf
from scipy import ndimage
import numpy as np
import os
import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import walnut_utils as wnut

meatextract = False
shellextract = True

src = '../clean/'
dsts = '../shells/'
dstm = '../meats/'

fs = 15
ssrange = [150,200,250,300]

Sx = [ np.s_[k,:,:] for k in ssrange ]
Sy = [ np.s_[:,k,:] for k in ssrange ]
Sz = [ np.s_[:,:,k] for k in ssrange ]

walnut_batch = sorted(glob.glob( src + '*/' ))

for bidx in range(len(walnut_batch)):
#for bidx in [4]:

    walnut_files = sorted(glob.glob(walnut_batch[bidx] + '*.tif'))
    bname = walnut_batch[bidx].split('/')[-2]
    wdsts = dsts + bname + '/'
    if not os.path.isdir(wdsts):
        os.makedirs(wdsts)
    wdstsd = wdsts + 'diagnostic/'
    if not os.path.isdir(wdstsd):
        os.makedirs(wdstsd)

    wdstm = dstm + bname + '/'
    if not os.path.isdir(wdstm):
        os.makedirs(wdstm)
    wdstmd = wdstm + 'diagnostic/'
    if not os.path.isdir(wdstmd):
        os.makedirs(wdstmd)

    for widx in range(len(walnut_files)):

        print(walnut_files[widx])
        img = tf.imread(walnut_files[widx])
        struc1 = ndimage.generate_binary_structure(img.ndim, 1)
        struc2 = ndimage.generate_binary_structure(img.ndim, 2)

        pa, fi = os.path.split(walnut_files[widx])
        fname = os.path.splitext(fi)[0]

        ################################################
        # Meat
        ################################################

        if meatextract:

            mimg = img.copy()

            eimg = ndimage.binary_erosion(mimg, struc2, 12, border_value=1)
            mimg = eimg*mimg

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

            dimg = ndimage.binary_dilation(eimg, struc2, 8, border_value=0)
            meat = dimg*img
            meat[meat < 120] = 0

            ################################################

            bmeat = meat.copy()
            bmeat[bmeat > 0] = 1
            bmeat = ndimage.binary_dilation(bmeat, struc2, 1, border_value=1).astype(meat.dtype)
            bmeat = ndimage.binary_erosion(bmeat, struc2, 1, border_value=0).astype(meat.dtype)

            fmeat = bmeat.copy()
            #fmeat = wnut.fill_component(fmeat)
            fmeat = ndimage.binary_fill_holes(fmeat, struc1)

            diffm = fmeat - bmeat
            labels,hist,numgeq1 = wnut.label_and_rearrange(diffm)

            if not numgeq1:
                labels, num = ndimage.label(diffm, structure=struc1)
                numholes = 0

            else:
                thr = 1e-4
                numholes = np.sum(np.sort(hist)[::-1]/np.sum(bmeat) > thr)
                print(numholes)

            bmeat[ labels > numholes ] = 1
            meat = bmeat*img
            print(hist)
            print(hist.max())
            print(np.sum(hist))

            ################################################

            fig, ax = plt.subplots(3,4,figsize=(16,14), sharex=True, sharey=True)
            ax = np.atleast_1d(ax).flatten()

            for i in range(4):
                ax[i].imshow(img[Sx[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)
            for i in range(4,8):
                ax[i].imshow(meat[Sx[i%4]], cmap='inferno', origin='lower', vmax=1, vmin=0)
            for i in range(8,len(ax)):
                ax[i].imshow(meat[Sx[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)

            fig.suptitle(bname + '/' + fname + ' : x-plane meat', fontsize=fs)
            fig.tight_layout()

            filename = wdstmd + bname + '_' + fname + '_proc_meat_x.jpg'
            plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
            plt.close()

            ################################################

            fig, ax = plt.subplots(3,4,figsize=(16,14), sharex=True, sharey=True)
            ax = np.atleast_1d(ax).flatten()

            for i in range(4):
                ax[i].imshow(img[Sy[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)
            for i in range(4,8):
                ax[i].imshow(meat[Sy[i%4]], cmap='inferno', origin='lower', vmax=1, vmin=0)
            for i in range(8,len(ax)):
                ax[i].imshow(meat[Sy[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)

            fig.suptitle(bname + '/' + fname + ' : y-plane meat', fontsize=fs)
            fig.tight_layout()

            filename = wdstmd + bname + '_' + fname + '_proc_meat_y.jpg'
            plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
            plt.close()

            ################################################

            filename = wdstm + bname + '_' + fname + '_meat.tif'
            tf.imwrite(filename, meat, photometric='minisblack', compress=3)

            snaps = wnut.collapse_dimensions(meat)
            wnut.plot_collapse_dimensions(snaps, bname+'_'+fname, 'meat', dst=wdstmd, writefig=True, display=False)

        ################################################
        # Shell
        ################################################

        if shellextract:

            filename = wdstm + bname + '_' + fname + '_meat.tif'
            meat = tf.imread(filename)

            timg = img.copy()
            timg[timg < 120] = 0
            timg = wnut.get_largest_element(timg)
            timg[timg > 0] = 1

            bmeat = meat.copy()
            bmeat[bmeat > 0] = 1

            bmeat = ndimage.binary_dilation(bmeat, struc2, 5, border_value=0)
            wmask = bmeat | timg

            wmask = ndimage.binary_dilation(wmask, struc2, 10, border_value=0)
            wmask = ndimage.binary_fill_holes(wmask)
            wmask = ndimage.binary_erosion(wmask, struc2, 5, border_value=1)

            fimg = img*wmask

            fimg[fimg > 175] = 0
            fimg[fimg == 0] = 2**(8*img.itemsize)-1

            nonmeat = ~meat.astype(np.bool)
            nonmeat = ndimage.binary_dilation(nonmeat, struc1, 2, border_value=0)
            simg = fimg*nonmeat
            simg[simg < 170] = 0

            uimg = ndimage.uniform_filter(simg, (3,3,3))
            uimg[uimg < 150] = 0

            eimg = ndimage.grey_erosion(uimg, (3,3,3), mode='constant', cval = 255)
            mcomp = wnut.get_largest_element(eimg)

            cshell = ndimage.binary_dilation(mcomp, struc2, 5, border_value=0)

            shell = img*cshell
            shell[shell < 150] = 0
            shell = wnut.get_largest_element(shell)

            ################################################

            bshell = shell.copy()
            bshell[bshell > 0] = 1
            bshell = ndimage.binary_dilation(bshell, struc2, 1, border_value=1).astype(shell.dtype)
            bshell = ndimage.binary_erosion(bshell, struc2, 1, border_value=0).astype(shell.dtype)

            fshell = bshell.copy()
            fshell = ndimage.binary_fill_holes(fshell, struc1)

            diffs = fshell - bshell
            sdiff = np.sum(diffs)
            sbshell = np.sum(bshell)
            print(sdiff, sbshell, sdiff < sbshell, sep='\t')

            thr = 1e-1

            if np.all(diffs == 0):
                shell = img*bshell
            elif sdiff < sbshell:
                shell = img*fshell
            else:
                labels,hist,numgeq1 = wnut.label_and_rearrange(diffs)
                holes = np.sum(np.sort(hist)[::-1]/np.sum(bmeat) > thr)
                print('holes', holes)
                print('max', np.max(hist))
                bshell[ labels > holes ] = 1
                shell = bshell*img

            ################################################

            fig, ax = plt.subplots(3,4,figsize=(16,14), sharex=True, sharey=True)
            ax = np.atleast_1d(ax).flatten()

            for i in range(4):
                ax[i].imshow(img[Sx[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)
            for i in range(4,8):
                ax[i].imshow(shell[Sx[i%4]], cmap='inferno', origin='lower', vmax=1, vmin=0)
            for i in range(8,len(ax)):
                ax[i].imshow(shell[Sx[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)

            fig.suptitle(bname + '/' + fname + ' : x-plane shell', fontsize=fs)
            fig.tight_layout()

            filename = wdstsd + bname + '_' + fname + '_proc_shell_x.jpg'
            plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
            plt.close()

            ################################################

            fig, ax = plt.subplots(3,4,figsize=(16,14), sharex=True, sharey=True)
            ax = np.atleast_1d(ax).flatten()

            for i in range(4):
                ax[i].imshow(img[Sy[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)
            for i in range(4,8):
                ax[i].imshow(shell[Sy[i%4]], cmap='inferno', origin='lower', vmax=1, vmin=0)
            for i in range(8,len(ax)):
                ax[i].imshow(shell[Sy[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)

            fig.suptitle(bname + '/' + fname + ' : y-plane shell', fontsize=fs)
            fig.tight_layout()

            filename = wdstsd + bname + '_' + fname + '_proc_shell_y.jpg'
            plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
            plt.close()

            ################################################

            fig, ax = plt.subplots(3,4,figsize=(16,14), sharex=True, sharey=True)
            ax = np.atleast_1d(ax).flatten()

            for i in range(4):
                ax[i].imshow(img[Sz[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)
            for i in range(4,8):
                ax[i].imshow(shell[Sz[i%4]], cmap='inferno', origin='lower', vmax=1, vmin=0)
            for i in range(8,len(ax)):
                ax[i].imshow(shell[Sz[i%4]], cmap='inferno', origin='lower', vmax=255, vmin=0)

            fig.suptitle(bname + '/' + fname + ' : z-plane shell', fontsize=fs)
            fig.tight_layout()

            filename = wdstsd + bname + '_' + fname + '_proc_shell_z.jpg'
            plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
            plt.close()

            ################################################

            filename = wdsts + bname + '_' + fname + '_shell.tif'
            tf.imwrite(filename, shell, photometric='minisblack', compress=3)

            snaps = wnut.collapse_dimensions(shell)
            wnut.plot_collapse_dimensions(snaps, bname+'_'+fname, 'shell', dst=wdstsd, writefig=True, display=False)

