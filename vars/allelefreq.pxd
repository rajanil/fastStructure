
import numpy as np
cimport numpy as np
cimport admixprop as ap
from cpython cimport bool

cdef class AlleleFreq:

    """
    Allele frequencies for SNPs and their relevant methods
    fall under this class. The simple prior over allele frequencies is set
    to be a "flat" Beta distribution. If the choice of prior is ``logistic'',
    the allele frequencies are assumed to be the logistic transform of
    a normally-distributed random variable with loci-specific means
    and population-specific variances.

    Arguments
        L : int
            number of loci
        K : int
            number of populations
        prior : str
            {'simple','logistic'}

    """

    cdef long L,K
    cdef double mintol
    cdef np.ndarray beta, gamma, var_beta, var_gamma, zetabeta, zetagamma, piA, F, mu, Lambda
    cdef list oldvar_beta, oldvar_gamma
    cdef str prior 

    cdef copy(self)
    cdef require(self)
    cdef _update_simple(self, np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi)
    cdef _update_logistic(self, np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi)
    cdef _unconstrained_solver(self, np.ndarray[np.float64_t, ndim=2] Dvarbeta, np.ndarray[np.float64_t, ndim=2] Dvargamma)
    cdef update(self, np.ndarray[np.uint8_t, ndim=2] G, AlleleFreq psi)
    cdef square_update(self, np.ndarray[np.uint8_t, ndim=2] G, AlleleFreq psi)
    cdef update_hyperparam(self, bool nolambda)
