import tifffile as tf
from scipy import ndimage, signal, optimize
import numpy as np
import os
import glob

import warnings
warnings.filterwarnings( "ignore")
from matplotlib import pyplot as plt

import unionfind as UF

def normalize_density(img, npz):
    resol = 2**(img.dtype.itemsize*8)
    rescale = np.round(npz(np.arange(resol, dtype=img.dtype)))
    rescale[rescale < 0] = 0
    rescale[rescale > resol - 1] = resol - 1
    rescale.astype(img.dtype)

    with np.nditer(img, flags=['external_loop'], op_flags=['readwrite']) as it:
        for x in it:
            x[...] = rescale[x]

    return img

def get_largest_element(comp, thr=0.1, outlabels=False):
    tot = np.sum(comp > 0)
    labels,num = ndimage.label(comp, structure=ndimage.generate_binary_structure(comp.ndim, 1))
    hist,bins = np.histogram(labels, bins=num, range=(1,num+1))
    argsort_hist = np.argsort(hist)[::-1]

    where = np.where(hist/tot > thr)[0] + 1
    print(num,'components\t',len(where),'preserved')
    print(np.sort(hist)[::-1][:20])

    mask = labels == where[0]
    for w in where[1:]:
        mask = mask | (labels == w)
    box0 = comp.copy()
    box0[~mask] = 0

    if outlabels:
        return box0, labels, where

    return box0

def label_and_rearrange(obinglands):
    numgeq1 = True

    struc = ndimage.generate_binary_structure(obinglands.ndim, 1)

    labels,num = ndimage.label(obinglands, structure=struc)
    print(num,'components')
    if num > 1:
        hist,bins = np.histogram(labels, bins=num, range=(1,num+1))
        argsort_hist = np.argsort(hist)[::-1]
        npz = np.hstack(([0],1+np.argsort(argsort_hist)))

        with np.nditer(labels, flags=['external_loop'], op_flags=['readwrite']) as it:
            for x in it:
                x[...] = npz[x]
    else:
        numgeq1 = False
        hist = np.arange(1)

    return labels, hist, numgeq1

def fill_component(comp, x=True, y=True, z=True):
    rcomp = comp.copy()
    rcomp[rcomp > 0] = 1

    if x:
        for k in range(rcomp.shape[0]):
            rcomp[k,:,:] = ndimage.binary_fill_holes(rcomp[k,:,:])
        print('Closed X')
    if y:
        for k in range(rcomp.shape[1]):
            rcomp[:,k,:] = ndimage.binary_fill_holes(rcomp[:,k,:])
        print('Closed Y')
    if z:
        for k in range(rcomp.shape[2]):
            rcomp[:,:,k] = ndimage.binary_fill_holes(rcomp[:,:,k])
        print('Closed Z')

    return rcomp

def clean_zeroes(img, pad=2):
    dim = img.ndim
    orig_size = img.size

    cero = np.arange(2*dim)

    for k in range(dim):
        ceros = np.all(img == 0, axis = (k, (k+1)%dim))

        for i in range(len(ceros)):
            if(~ceros[i]):
                break
        for j in range(len(ceros)-1, 0, -1):
            if(~ceros[j]):
                break
        cero[k] = i
        cero[k+dim] = j+1
    for i in range(dim):
        cero[i] -= 2
    for i in range(dim, len(cero)):
        cero[i] += 2
    cero[cero < 0] = 0
    img = img[cero[1]:cero[4], cero[2]:cero[5], cero[0]:cero[3]]

    print(round(100-100*img.size/orig_size),'% reduction from input')

    return img, cero

def collapse_dimensions(img):
    snaps = []
    for i in range(img.ndim):
        snaps.append(np.sum(img, axis=i))
    return snaps

def collapse_dimensions_max(img):
    snaps = []
    for i in range(img.ndim):
        snaps.append(np.max(img, axis=i))
    return snaps

def plot_collapse_dimensions(snaps, bname='bname', tissue='tissue', display=False, writefig=False, dst='./'):
    fig, ax = plt.subplots(1,len(snaps),figsize=(4*len(snaps),5))
    for i in range(len(snaps)):
        ax[i].imshow(snaps[i], cmap='inferno', origin='lower');
    plt.suptitle(bname + ' ' + tissue + ' collapse', fontsize=20);
    plt.tight_layout()

    if writefig:
        filename = dst + bname + '_' + '_'.join(tissue.split(' ')) + '.jpg'
        plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight');
    if not display:
        plt.close();
    return fig,ax


########################################################################
########################################################################

def ellipsoid(lamb, phi, ax=3, ay=2, az=1):

    ex2 = (ax*ax - az*az)/(ax*ax)
    ey2 = (ay*ay - az*az)/(ay*ay)
    ee2 = (ax*ax - ay*ay)/(ax*ax)

    nu = ax/(np.sqrt(1 - ex2*np.sin(phi)**2 - ee2*(np.cos(phi)*np.sin(lamb))**2))

    return np.array([nu*np.cos(phi)*np.cos(lamb),
                     nu*(1-ee2)*np.cos(phi)*np.sin(lamb),
                     nu*(1-ex2)*np.sin(phi)])

def ell_algebraic_fit_3d(datapoints, k=4):

    x,y,z = datapoints
    C1 = np.diag([-1, -1, -1, -k, -k, -k])
    C1[0,1:3] = .5*k - 1
    C1[1:3, 0] = .5*k - 1
    C1[1,2] = .5*k - 1
    C1[2,1] = .5*k - 1
    D = np.vstack((x**2, y**2, z**2, 2*x*y, 2*x*z, 2*y*z, 2*x, 2*y, 2*z, np.ones_like(x)))

    DDT = D@D.T
    S11 = DDT[:6, :6]
    S12 = DDT[:6, 6:]
    S22 = DDT[6:, 6:]

    M = (np.linalg.inv(C1))@(S11 - S12@np.linalg.inv(S22)@(S12.T))
    eigvals, eigvecs = np.linalg.eig(M)
    v1 = eigvecs[:, np.argmax(eigvals)]
    v2 = -np.linalg.inv(S22)@(S12.T)@v1

    v = -np.hstack((v1,v2))/v2[-1]

    a,b,c,f,g,h, p,q,r, d = v
    I = a + b + c
    J = a*b + b*c + c*a - f*f/4 - g*g/4 - h*h/4
    aux = k*J - I*I

    return v, aux > 0

def get_ell3d_params_from_vector(c, ax_guess=None):
    if ax_guess is None:
        ax_guess = [0,1,2]
    ax_guess = np.asarray(ax_guess)

    Cmat = np.array([
        [2*c[0], c[3]  , c[4]  ],
        [c[3],   2*c[1], c[5]  ],
        [c[4],   c[5]  , 2*c[2]]
    ])
    t = -np.linalg.inv(Cmat) @ c[6:9]

    D = 1 + np.sum(c[:3]*(t**2)) + c[3]*t[0]*t[1] + c[4]*t[0]*t[2] + c[5]*t[1]*t[2]

    Q = D*np.linalg.inv(.5*Cmat)
    eigvals, R = np.linalg.eigh(Q)
    axis = np.sqrt(eigvals)[ax_guess]
    for i,cc in enumerate(np.column_stack((np.argmax(np.abs(R), axis = 0), [0,1,2]))):
        if(R[tuple(cc)] < 0 ):
            R[:,i] = -R[:,i]
    R = R[:, ax_guess]
    R = R.T
    t = t

    thetax = np.arctan(-R[2,1]/R[2,2])
    thetay = np.arctan( R[2,0]/(np.sqrt(R[0,0]**2 + R[1,0]**2)))
    thetaz = np.arctan(-R[1,0]/R[0,0])

    ell_params = {'center': t,
                  'axes': axis,
                  'rotation': R,
                  'theta': np.array([thetax, thetay, thetaz])}
    return ell_params

########################################################################
########################################################################

def find_tip_min(coords, x=0, y=1, z=2):
    maxes = np.min(coords, axis = 1)
    max_vox = coords[:, coords[z] == maxes[z]]
    if len(max_vox) > 1 :
        max_vox = np.mean(max_vox, axis=1)
        #maxesz = np.min(max_vox, axis = 1)
        #max_vox = max_vox[:, max_vox[y] == maxesz[y]]

        #if len(max_vox) > 1:
        #    maxesy = np.min(max_vox, axis = 1)
        #    max_vox = max_vox[:, max_vox[x] == maxesy[x]]

    return np.squeeze(max_vox)

def find_tip_max(coords, x=0, y=1, z=2):
    maxes = np.max(coords, axis = 1)
    max_vox = coords[:, coords[z] == maxes[z]]
    if len(max_vox) > 1 :
        max_vox = np.mean(max_vox, axis=1)
        #maxesz = np.max(max_vox, axis = 1)
        #max_vox = max_vox[:, max_vox[y] == maxesz[y]]

        #if len(max_vox) > 1:
        #    maxesy = np.max(max_vox, axis = 1)
        #    max_vox = max_vox[:, max_vox[x] == maxesy[x]]

    return np.squeeze(max_vox)

def plot_3Dprojections(seed, title='title', markersize=2, alpha=1, mk = '.', origin=[0,0,0], writefig=False, dst='./', dpi=125, display=True):
    axes = ['X','Y','Z']
    fig, ax = plt.subplots(1,3,figsize=(15,4))

    for i in range(3):
        proj = []
        for j in range(3):
            if j != i:
                proj.append(j)
        ax[i].scatter(seed[proj[1]], seed[proj[0]], s=markersize, color='y', alpha=alpha, marker=mk)
        ax[i].set_ylabel(axes[proj[0]])
        ax[i].set_xlabel(axes[proj[1]])
        ax[i].set_title(axes[i] + ' Projection')
        ax[i].set_aspect('equal', 'datalim');

        ax[i].axvline(origin[proj[1]], c='r',alpha=0.2)
        ax[i].axhline(origin[proj[0]],c='r',alpha=0.2)

    fig.suptitle(title, y=0.95, fontsize=20)
    plt.tight_layout();

    if writefig:
        filename = '_'.join(title.split(' '))
        plt.savefig(dst + filename + '.png', dpi=dpi, format='png', bbox_inches='tight',
                    facecolor='white', transparent=False)
    if not display:
        plt.close();

    return fig, ax

########################################################################
########################################################################

def ell_algebraic_fit_2d(X,Y, verbose=False):
    A = np.column_stack([X**2, X * Y, Y**2, X, Y])
    b = np.ones_like(X)
    x = np.linalg.lstsq(A, b)[0].squeeze()
    if verbose:
        print('The ellipse is given by {0:.3}x^2 + {1:.3}xy+{2:.3}y^2+{3:.3}x+{4:.3}y = 1'.format(*x))

    return x

def get_ell2d_params_from_vector(V):
    A,B,C,D,E = V
    F = -1

    det = B*B - 4*A*C
    foo = np.sqrt((A-C)**2 + B**2)
    bar = -np.sqrt(2*(A*E*E + C*D*D + B*D*E + det*F))/det

    a = bar*np.sqrt(A + C + foo)
    b = bar*np.sqrt(A + C - foo)

    origx = (2*C*D - B*E)/det
    origy = (2*A*E - B*D)/det

    if B != 0:
        theta = np.arctan((C-A-foo)/B)
    elif A < C:
        theta = 0
    else:
        theta = np.pi*.5

    rotation = np.array([[np.cos(theta), -np.sin(theta)],
                         [np.sin(theta),  np.cos(theta)]])

    return {'axes':np.array([a,b]), 'origin': np.array([origx, origy]), 'theta':theta, 'rot': rotation}

def Ft(t, y0=1, y1=1, e0=5, e1=1):
    return (e0*y0/(t + e0*e0))**2 + (e1*y1/(t + e1*e1))**2 - 1

def d_ell0_point(y0=1, y1=1, e0=5, e1=1):
    x0,x1 = 0,0
    distance = -1
    sgn = 0

    if y1 > 0:
        if y0 > 0:
            t0 = -e1*e1 + e1*y1
            t1 = -e1*e1 + np.sqrt(e0*e0*y0*y0 + e1*e1*y1*y1)

            bracket = [t0, t1]
            sol = optimize.root_scalar(Ft, args=(y0,y1,e0,e1), method='bisect', bracket=bracket, maxiter=10000)

            if sol.converged:
                x0 = e0*e0*y0/(sol.root + e0*e0)
                x1 = e1*e1*y1/(sol.root + e1*e1)
                distance = np.sqrt((x0-y0)**2 + (x1-y1)**2)
                sgn = np.sign(sol.root)

        else: # y0 == 0
            x0 = 0
            x1 = e1
            distance = np.abs(y1-e1)
            sgn = np.sign(y1-e1)

    else: # y1 == 0
        if y0 < (e0*e0 - e1*e1)/e0:
            x0 = e0*e0*y0/(e0*e0 - e1*e1)
            x1 = e1*np.sqrt(1 - x0*x0/(e0*e0))
            distance = np.sqrt((x0-y0)**2 + x1*x1)
            sgn = np.sign(x0-y0)

        else:
            x0 = e0
            x1 = 0
            distance = np.abs(y0 - e0)
            sgn = np.sign(y0 - e0)

    return distance, sgn, x0, x1

def d_ell_point(y0,y1, pdict):
    #y0,y1 = pdict['rot'] @ (Q - pdict['origin'])
    #y0, y1 = Q

    if (y0 >= 0) and (y1 >= 0):
        d, s, x0, x1 = d_ell0_point(y0,y1, *pdict['axes'])
        return d, s, x0, x1
    elif (y0 >= 0) and (y1 < 0):
        d, s, x0, x1 = d_ell0_point(y0, -y1, *pdict['axes'])
        return d, s, x0, -x1
    elif (y0 < 0) and (y1 < 0):
        d, s, x0, x1 = d_ell0_point(-y0, -y1, *pdict['axes'])
        return d, s, -x0, -x1
    elif (y0 < 0) and (y1 >= 0):
        d, s, x0, x1 = d_ell0_point(-y0, y1, *pdict['axes'])
        return d, s, -x0, x1

def even_space_ell0(N=100, a=2, b=1):
    # Following https://math.stackexchange.com/questions/2093569/points-on-an-ellipse

    e = 1 - b*b/(a*a)
    t = np.linspace(0, 2*np.pi, N)
    theta = t + (e**2/8 + e**4/16 + 71*e**6/2048)*np.sin(2*t) + (5*e**4 + 5*e**6)*np.sin(4*t)/256 + 29*e**6*np.sin(6*t)/6144

    x = a*np.cos(theta)
    y = b*np.sin(theta)

    return np.vstack((theta, x, y))
