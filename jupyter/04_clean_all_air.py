import tifffile as tf
from scipy import ndimage
import numpy as np

import os
import glob


src = '../clean2/'
dstc = '../clean2/'

struc1 = ndimage.generate_binary_structure(3, 1)
struc2 = ndimage.generate_binary_structure(3, 2)
struc2 = ndimage.generate_binary_structure(3, 3)

walnut_batch = sorted(glob.glob( src + '*/' ))

for bidx in range(len(walnut_batch)):

    walnut_files = sorted(glob.glob(walnut_batch[bidx] + '*.tif'))
    bname = walnut_batch[bidx].split('/')[-2]

    wdstc = dstc + bname + '/'
    if not os.path.isdir(wdstc):
        os.makedirs(wdstc)


    for widx in range(len(walnut_files)):
        print(walnut_files[widx])
        img = tf.imread(walnut_files[widx])

        pa, fi = os.path.split(walnut_files[widx])
        fname = os.path.splitext(fi)[0]

        sz = 7
        dist1 = ndimage.distance_transform_cdt(img, metric='taxicab')

        mask1 = (img > 130) | (dist1 > 10)
        mask1 = ndimage.binary_dilation(mask1, struc1, 7, border_value=0)
        mask1 = ndimage.binary_fill_holes(mask1)
        mask1 = ndimage.binary_erosion(mask1, struc1, 7, border_value=1)
        clean = mask1*img


        filename = wdstc + bname + '_' + fname + '.tif'
        filename = walnut_files[widx]
        print(filename)
        tf.imwrite(filename, clean, photometric='minisblack', compress=3)
