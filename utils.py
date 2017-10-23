#!/usr/bin/env python3
# © 2014 Nicolai Waniek <rochus at rochus dot net>, see LICENSE file for details

"""Utility functions."""

from numpy import mean, std, ones, conj, prod, isnan, sqrt, dot, copy, \
        linspace, zeros, amax, multiply, sum, mod, floor, array
from scipy.signal import correlate2d
from scipy.ndimage.filters import convolve, gaussian_filter
from scipy.ndimage.interpolation import rotate


def sumsq(x):
    return sum(x * conj(x))


def toXY(idx, N):
    """return an array([Row,Col]) of the coordinate describing the index Idx
    into an NxN matrix"""

    return array([floor(idx / N), mod(idx, N)], dtype=int)


def normcorr2d(A, B=None):
    """normalized correlation of a 2D input array A with mask B. If B is
    not specified, the function computes the normalized utocorrelation
    of A."""

    # TODO: check the dimensions

    if B == None:
        B = A

    N = A.shape[0] * A.shape[1]
    R = A - mean(A)
    S = B - mean(B)
    T = correlate2d(R, S, mode='full') / (N * std(A) * std(B))
    return T


def gridness_score(R):
    """Compute the Gridness Score of a rate map R"""

    dim0 = R.shape[0]
    print(dim0)
    cntr = dim0 // 2

    # create a ring filter to search for the six closest fields
    in_ra = 12;
    out_ra = 30;
    RingFilt = zeros((dim0, dim0))
    for i in range(dim0):
        for j in range(dim0):
                cntr_i = (cntr - i)**2
                cntr_j = (cntr - j)**2
                dist = cntr_i + cntr_j
                if in_ra**2 <= dist and dist <= out_ra**2:
                    RingFilt[i, j] = 1

    Dong = multiply(R, RingFilt)
    Dong = Dong[cntr - out_ra - 1:cntr + out_ra + 1, cntr - out_ra - 1 : cntr + out_ra + 1]
    nx, ny = Dong.shape[0], Dong.shape[1]

    corrot = zeros(180)
    for idrot in range(180):
        rot = rotate(Dong, idrot, reshape=False)
        A = Dong
        B = rot

        dotAA = sum(A*A, axis=0)
        dotA0 = sum(A*ones((nx, ny)), axis=0)
        dotBB = sum(B*B, axis=0)
        dotB0 = sum(B*ones((nx, ny)), axis=0)
        dotAB = sum(A*B, axis=0)

        corrot[idrot] = (nx * ny * sum(dotAB) - sum(dotA0) * sum(dotB0)) / \
                 (sqrt(nx * ny * sum(dotAA) - sum(dotA0)**2)*sqrt(nx*ny*sum(dotBB)) - sum(dotB0)**2)

    #plt.imshow(Dong, interpolation='None')
    #plt.plot(corrot)

    gridscore = min(corrot[59],corrot[119]) - max(corrot[29], corrot[89], corrot[149])
    return gridscore


# TODO: develop a version where array may be multidimensional, and idx
# is thus a tuple containing the closest value for each broadcast
# dimension
def find_nearest(arr, val):
    """Find the entry in arr which is closes to val"""
    idx = (abs(arr - val)).argmin()
    return idx, arr[idx]



def spikes_to_ratemap(spikes, minxy=-0.5, maxxy=0.5, resolution=50, filter=True):
    """Convert an array spikes of the form Mx2, with M rows of spikes
    giving x/y coordinates into a rate map that is normalized to be in
    the range [0,1]. The additional option to filter will filter the resulting map with a gaussian and sigma=2.0"""

    mapping = linspace(minxy, maxxy, resolution)
    R = zeros((resolution, resolution))

    for spike in spikes:
        x,_ = find_nearest(mapping, spike[0])
        y,_ = find_nearest(mapping, spike[1])
        R[y,x] += 1
    R[isnan(R)] = 0

    if filter:
        R = gaussian_filter(R, (2.0, 2.0))
    R = R / amax(R)
    return R


# Testing Functions

# TODO: other place? something...
def test_autocorr():
    from numpy.random import random, randn
    from numpy.testing import assert_equal

    a = 10 * randn(100, 100)
    auto = normxcorr(a)
    add_in = normxcorr(a, -a)
    assert_equal(auto, -add_in)
