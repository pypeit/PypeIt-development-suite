import numpy as np
from pypeit.traceslits import TraceSlits
from pypeit.spectrographs.util import load_spectrograph
from pypeit import ginga
from scipy import interpolate
from scipy import ndimage
from pypeit import msgs
from skimage.transform import resize


# if not (np.mod(a.shape[0], newshape[0]) == 0) or (np.mod(newshape[0], a.shape[0]) == 0):
#    msgs.error('a.shape[0] and newshape[0] are not interger multiples of each other')

# if not (np.mod(a.shape[1], newshape[1]) == 0) or (np.mod(newshape[1], a.shape[1]) == 0):
#    msgs.error('a.shape[1] and newshape[1] are not interger multiples of each other')


def rebin(a, newshape):
    '''Rebin an array to a new shape using slicing. This routine is taken from:
    https://scipy-cookbook.readthedocs.io/items/Rebinning.html
    '''

    if not len(a.shape) == len(newshape):
        msgs.error('Dimension of a image does not match dimension of new requested image shape')

    slices = [slice(0, old, float(old) / new) for old, new in zip(a.shape, newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')  # choose the biggest smaller integer index
    return a[tuple(indices)]


def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array(a.shape)
    ndims = len(a.shape)
    if len(newdims) != ndims:
        print("[congrid] dimensions error. This routine currently only support rebinning to the same number of dimensions.")
        return None
    newdims = np.asarray(newdims, dtype=float)
    dimlist = []

    if method == 'neighbour':
        for i in range(ndims):
            base = np.indices(newdims.astype(int))[i]
            dimlist.append((old[i] - m1) / (newdims[i] - m1) \
                           * (base + ofs) - ofs)
        cd = np.array(dimlist).round().astype(int)
        from IPython import embed
        embed()
        newa = a[list(cd)]
        return newa

    elif method in ['nearest', 'linear']:
        # calculate new dims
        for i in range(ndims):
            base = np.arange(newdims[i])
            dimlist.append((old[i] - m1) / (newdims[i] - m1) \
                           * (base + ofs) - ofs)
        # specify old dims
        olddims = [np.arange(i, dtype=np.float) for i in list(a.shape)]

        # first interpolation - for ndims = any
        mint = interpolate.interp1d(olddims[-1], a, kind=method)
        newa = mint(dimlist[-1])

        trorder = [ndims - 1] + range(ndims - 1)
        for i in range(ndims - 2, -1, -1):
            newa = newa.transpose(trorder)

            mint = interpolate.interp1d(olddims[i], newa, kind=method)
            newa = mint(dimlist[i])

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose(trorder)

        return newa
    elif method in ['spline']:
        oslices = [slice(0, j) for j in old]
        oldcoords = np.ogrid[oslices]
        nslices = [slice(0, j) for j in list(newdims)]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.rank(newcoords))
        # make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type. Currently only \'neighbour\', \'nearest\',\'linear\',", \
        "and \'spline\' are supported.")
        return None


TSlits = TraceSlits(None, None)
masterfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT_old/Keck_NIRES/NIRES/MF_keck_nires/MasterTrace_A_15_01'
tset_slits = TSlits.load_master(masterfile)
spectrograph = load_spectrograph('keck_nires')
slitmask_orig = spectrograph.slitmask(tset_slits)
slitmask_orig = slitmask_orig[1:,1:]
#slitmask_new = (np.round(congrid(slitmask_orig.astype(np.float64), (2048,2048), method='neighbour'))).astype(slitmask_orig.dtype)
newshape = (2047*2,1023)
slitmask_new = rebin(slitmask_orig, newshape)
slitmask_old = rebin(slitmask_new,slitmask_orig.shape)
slitmask_new1 = ((np.round(resize(slitmask_orig.astype(np.integer), newshape, preserve_range=True, order=0))).astype(np.integer)).astype(slitmask_orig.dtype)
slitmask_old1 = ((np.round(resize(slitmask_new1.astype(np.integer), slitmask_orig.shape, preserve_range=True, order=0))).astype(np.integer)).astype(slitmask_orig.dtype)
