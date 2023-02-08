import tifffile as tf
from scipy import ndimage

import numpy as np
import pandas as pd

import os
import glob
import argparse

import walnut_utils as wnut


def main():

    resol = 0.0759
    Tvals = [4,8,16,32,64]
    Tmax = max(Tvals)
    mborder = ndimage.generate_binary_structure(3, 3).astype(int)
    mborder[mborder == 1] = -1
    mborder[1,1,1] = -np.sum(mborder) - 1
    buckets = [None for i in range(4)]

    parser = argparse.ArgumentParser(description='Normalize density values of walnuts')

    parser.add_argument('src', metavar='raw_walnut_src', type=str, help='path to raw walnut images')
    parser.add_argument('dst', metavar='clean_img_dst', type=str, help='path to store clean images')
    parser.add_argument('bname', metavar='scan_id', type=str, help='walnut batch scan id')
    args = parser.parse_args()

    src = args.src
    dst = args.dst
    bname = args.bname

    # src = '../hpcc/'
    # dst = '../hpcc/topology/'

    csrc = src + 'clean/'
    rsrc = src + 'rotated/'
    tsrc = src + 'watershed/'

    walnut_files = sorted(glob.glob(csrc + bname + '/*.tif'))
    direction_files = sorted(glob.glob(dst + '*.csv'))
    Udirections = [np.loadtxt(direction_file, delimiter=',') for direction_file in direction_files]

    wdst = dst + bname + '/'
    if not os.path.isdir(wdst):
        os.makedirs(wdst)

    for widx in range(len(walnut_files)):
        pa, fi = os.path.split(walnut_files[widx])
        fname = os.path.splitext(fi)[0]
        shell_file = tsrc + bname + '/' + fname + '_shell.tif'
        kernel_file = tsrc + bname + '/' + fname + '_meat.tif'
        print(fname)

        for imgfile, imgname in zip([walnut_files[widx], shell_file, kernel_file], ['walnut', 'shell', 'kernel']):

            if len(direction_files) == len(glob.glob(wdst + fname + '_' + imgname + '_*_t{:02d}.csv'.format(Tmax))):
                print('All ECTs already computed for whole', imgname)
                print(fname)

            else:

                img = tf.imread(imgfile)
                img = np.pad(img, 1)

                bimg = img.copy().astype(int)
                bimg[bimg > 0]  = 1
                border = ndimage.convolve(bimg, mborder, mode='constant', cval=0)
                border[border < 0] = 0
                bimg = bimg.astype(float).ravel()

                orig_datapoints = np.asarray(np.nonzero(border))
                flat_datapoints = np.flatnonzero(border).astype(np.uint64)

                # # Complexify

                mask = np.zeros(border.size, dtype=np.uint64)
                mask[flat_datapoints] = flat_datapoints

                edges, faces, cubes = wnut.complexify(flat_datapoints, mask, img.shape)

                nv = len(flat_datapoints)
                ne = len(edges)
                nf = len(faces)
                nc = len(cubes)
                print('\n----\n----',imgname.upper(),'\n----')
                print('Verts:\t',nv)
                print('Edges:\t',ne)
                print('Faces:\t',nf)
                print('Cubes:\t',nv)
                print('\nEC:\t',nv-ne+nf-nc)

                # # Rotate

                filename = rsrc + bname + '/' + fname + '_rotation.csv'
                data = np.loadtxt(filename, delimiter=',')

                wmean = data[0]
                rotxyz = data[10:13]
                coords = resol*(rotxyz @ (orig_datapoints - wmean.reshape(-1,1)))

                for ddx in range(len(direction_files)):
                    udirections = Udirections[ddx]
                    filename = wdst + fname + '_' + imgname + '_d{:04d}_t{:02d}.csv'.format(len(udirections), Tmax)
                    print(filename)
                    if os.path.isfile(filename):
                        print(filename, 'already computed')

                    else:
                        # # ECT

                        T = Tmax

                        ect = np.zeros(T*len(udirections), dtype=int)

                        for i in range(len(udirections)):
                            direction = udirections[i]
                            heights = np.sum(coords*direction.reshape(-1,1), axis=0)
                            bimg[flat_datapoints] = heights
                            minh = np.min(heights)
                            maxh = np.max(heights)
                            buckets[0],_ = np.histogram(heights, bins = T, range=(minh, maxh))
                            buckets[1],_ = np.histogram(np.max(bimg[edges], axis=1), bins=T, range=(minh, maxh))
                            buckets[2],_ = np.histogram(np.max(bimg[faces], axis=1), bins=T, range=(minh, maxh))
                            buckets[3],_ = np.histogram(np.max(bimg[cubes], axis=1), bins=T, range=(minh, maxh))

                            ecc = buckets[0] - buckets[1] + buckets[2] - buckets[3]
                            ecc = np.cumsum(ecc)
                            ect[i*T : (i+1)*T] = ecc

                        foo = pd.DataFrame(ect).T
                        foo.to_csv(filename, header=False, index=False)

    return 0

if __name__ == '__main__':
    main()
