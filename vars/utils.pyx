import numpy as np
cimport numpy as np
import random

"""
Collection of useful functions.
"""

EPS = np.finfo(np.double).tiny

insum = lambda x,axes: np.apply_over_axes(np.sum,x,axes)
nplog = lambda x: np.nan_to_num(np.log(x))
OR = lambda x,y: np.logical_or(x,y)

def outsum(np.ndarray[np.float64_t, ndim=2] arr):

    """
    Fast summation over the 0-th axis.
    """

    cdef list shape
    cdef np.ndarray thesum
    thesum = sum([a for a in arr])
    thesum = sum([a for a in arr])
    thesum = thesum.reshape(1,thesum.size)
    return thesum

def random_combination(iterable, r):

    """
    Random selection from itertools.combinations(iterable, r)
    """

    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(xrange(n), r))
    return tuple(pool[i] for i in indices)

def random_partition(list shape, list nonlst, int n):

    cdef int num, d, j
    cdef double frac
    cdef tuple i
    cdef list index, indices, validindices, masks

    frac = 0.01
    num = min([1000, int(frac*(shape[0]*shape[1]-len(nonlst)))])
    indices = [[(np.random.randint(0,shape[0]), np.random.randint(0,shape[1])) \
        for d in xrange(2*num)] for j in xrange(n)]
    validindices = [list(set(index).difference(set(nonlst))) for index in indices]
    masks = [(np.array([i[0] for i in index[:num]], dtype='int'), \
        np.array([i[1] for i in index[:num]], dtype='int')) for index in validindices]

    return masks
