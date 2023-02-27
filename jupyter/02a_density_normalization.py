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

fs = 15
lw = 3
thr = 100

y = np.array([40.5,172.,223.5])
anchors = np.array([0.52268214, 0.92105655, 0.98927759])

parser = argparse.ArgumentParser(description='Normalize density values of walnuts')

parser.add_argument('src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')

parser.add_argument('dst', metavar='clean_img_dst', type=str, help='path to store clean images')

parser.add_argument('bname', metavar='scan_id', type=str, help='walnut batch scan id')

args = parser.parse_args()

# src = '../raw/
# dst = '../clean/'
# bname = '2014SBa_R1_T25'

dst = args.dst
src = args.src
bname = args.bname

walnut_files = sorted(glob.glob(src + bname + '/*.tif'))

print(src + bname + '/*.tif')
print(walnut_files)

wdst = dst + bname + '/'

if not os.path.isdir(wdst + 'normalization/'):
    os.makedirs(wdst + 'normalization/')

print(wdst)

for widx in range(len(walnut_files)):
    print(walnut_files[widx])

    img = tf.imread(walnut_files[widx])

    if np.max(img) > 256:
        img = img//256
        img = img.astype(np.uint8)

    pa, fi = os.path.split(walnut_files[widx])
    fname = os.path.splitext(fi)[0]
    hist0,bins = np.histogram(img, bins=2**(img.dtype.itemsize*8),range=(0,2**(img.dtype.itemsize*8)))
    fhist = ndimage.median_filter(hist0, size=5, mode='constant', cval=0)
    cumul = np.cumsum(fhist)
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

    npz = P.Polynomial.fit(x,y,1)
    np.savetxt(wdst + 'normalization/normalization_coefficients' + fname + '.csv', npz.convert().coef, delimiter = ',')

    ################################################

    aimg = img.copy()
    aimg = wnut.normalize_density(aimg, npz)

    struc = ndimage.generate_binary_structure(img.ndim, 2)
    cimg = aimg.copy()
    cimg[cimg < thr] = 0
    oimg = wnut.fill_component(cimg)
    oimg = wnut.get_largest_element(oimg)

    struc = ndimage.generate_binary_structure(img.ndim, 1)

    eimg = ndimage.binary_erosion(oimg, struc, 3, border_value=1)
    eimg = ndimage.binary_dilation(eimg, struc, 3, border_value=0)

    cimg = (eimg > 0)*aimg

    ################################################

    k = 250
    ss = np.s_[k,:, :]

    fig, ax = plt.subplots(2,2,figsize=(15,15), sharex=True, sharey=True)
    ax = np.atleast_1d(ax).flatten()

    i = 0
    ax[i].imshow(img[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
    ax[i].set_title('img', fontsize=fs)
    i = 1
    ax[i].imshow(aimg[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
    ax[i].set_title('normalized', fontsize=fs)
    i = 2
    ax[i].imshow(oimg[ss], cmap='inferno', origin='lower', vmax=1, vmin=0)
    ax[i].set_title('components filled', fontsize=fs)
    i = 3
    ax[i].imshow(cimg[ss], cmap='inferno', origin='lower', vmax=255, vmin=0)
    ax[i].set_title('air removed', fontsize=fs)

    fig.tight_layout()
    fig.suptitle(bname + '/' + fname, fontsize=fs)
    filename = wdst + 'normalization/air_removal_' + fname + '.jpg'
    plt.savefig(filename, dpi=100, format='jpg', bbox_inches = 'tight', pil_kwargs={'optimize':True})
    plt.close()

    ################################################

    walnut, cero = wnut.clean_zeroes(cimg)
    np.savetxt(wdst + 'normalization/clean_zeroes' + fname + '.csv', cero, fmt='%d', delimiter = ',')

    bwalnut = walnut.copy()
    bwalnut[bwalnut > 0] = 1

    fwalnut = walnut.copy()
    fwalnut[fwalnut > 0] = 1
    fwalnut = ndimage.binary_fill_holes(fwalnut)

    diff = fwalnut - bwalnut

    walnut[diff > 0] = 1

    dist1 = ndimage.distance_transform_cdt(walnut, metric='taxicab')
    mask1 = (walnut > 130) | (dist1 > 15)
    mask1 = ndimage.binary_dilation(mask1, struc, 7, border_value=0)
    mask1 = ndimage.binary_fill_holes(mask1)
    mask1 = ndimage.binary_erosion(mask1, struc, 7, border_value=1)
    clean = mask1*walnut

    filename = wdst + bname + '_' + fname + '.tif'
    tf.imwrite(filename, clean, photometric='minisblack', compress=3)

    snaps = wnut.collapse_dimensions(clean)
    wnut.plot_collapse_dimensions(snaps, bname+'_'+fname, 'whole walnut', dst=wdst+'normalization/', writefig=True)

    snaps = wnut.collapse_dimensions_max(clean)
    wnut.plot_collapse_dimensions(snaps, bname+'_'+fname, 'walnut shell', dst=wdst+'normalization/', writefig=True)
