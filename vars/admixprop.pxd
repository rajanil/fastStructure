import numpy as np
cimport numpy as np
cimport allelefreq as af

cdef class AdmixProp:

    """
    Admixture proportions for all samples and their relevant methods
    fall under this class. The prior over admixture proportions is set
    to be a symmetric Dirichlet distribution with parameter 1/K.

    Arguments
        N : int
            number of samples
        K : int
            number of populations

    """

    cdef long N,K
    cdef np.ndarray alpha, var, xi
    cdef list oldvar

    cdef copy(self)
    cdef require(self)
    cdef update(self, np.ndarray[np.uint8_t, ndim=2] G, af.AlleleFreq pi)
    cdef square_update(self, np.ndarray[np.uint8_t, ndim=2] G, af.AlleleFreq pi)
