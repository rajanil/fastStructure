import numpy as np
cimport numpy as np
from cpython cimport bool
import utils

ctypedef np.uint8_t uint8_t

cdef extern from "admixprop.h":
    void Q_update( uint8_t* G, double* zetabeta, double* zetagamma, double* xi, double* new_var, long N, long L, long K )

cdef class AdmixProp:
    """Admixture proportions in samples.

    Arguments
        N : int
            number of samples
        K : int
            number of populations

    """

    cdef long N,K
    cdef np.ndarray alpha, var, xi
    cdef list oldvar

    def __cinit__(self, long N, long K):

        self.N = N
        self.K = K

        # Initializing hyperparameters
        self.alpha = 1./K*np.ones((1,K))

        # Initializing variational parameters for admixture proportions
        self.var = np.ones((N,K)) + 0.1*np.random.rand(N,K)
#        self.var = np.array([np.arange(1,K+1) for n in xrange(N)])
        self.xi = np.exp(digamma(self.var)-digamma(insum(self.var,[1])))
        self.oldvar = []

    cdef copy(self):

        cdef AdmixProp newinstance
        newinstance = AdmixProp(self.N, self.K)
        newinstance.var = self.var.copy()
        newinstance.xi = self.xi.copy()
        newinstance.oldvar = []

        return newinstance

    cdef require(self):

        self.var = np.require(self.var, dtype=np.float64, requirements='C')
        self.xi = np.require(self.xi, dtype=np.float64, requirements='C')

    cdef update(self, np.ndarray[np.uint8_t, ndim=2] G, Pi pi):
        """Update variational parameters of distribution over admixture proportions in the E step.

        Arguments
            G : numpy array

            pi : instance of `Pi`

        """

        self.var = np.zeros((self.N,self.K), dtype=np.float64)
        psi_update(<np.uint8_t*>G.data, <double*> pi.zetabeta.data, <double*> pi.zetagamma.data, <double*> self.xi.data, <double*> self.var.data, self.N, pi. L, self.K)

        if np.isnan(self.var).any():
            self.var = self.oldvar[-1]
        else:
            self.var = self.alpha + self.var
            self.xi = np.exp(digamma(self.var)-digamma(insum(self.var,[1])))
        self.require()

    cdef square_update(self, np.ndarray[np.uint8_t, ndim=2] G, Pi pi):

        cdef long step
        cdef bool a_ok
        cdef np.ndarray R, V
        cdef double a

        self.oldvar = [self.var.copy()]
        for step from 0 <= step < 2:
            self.update(G, pi)
            self.oldvar.append(self.var.copy())

        R = self.oldvar[1] - self.oldvar[0]
        V = self.oldvar[2] - self.oldvar[1] - R
        a = -1.*np.sqrt((R*R).sum()/(V*V).sum())

        if a>-1:
            a = -1.

        a_ok = False
        while not a_ok:
            self.var = (1+a)**2*self.oldvar[0] - 2*a*(1+a)*self.oldvar[1] + a**2*self.oldvar[2]
            if (self.var<=0).any():
                a = (a-1)/2.
                if np.abs(a+1)<1e-4:
                    a = -1.
            else:
                a_ok = True

        # if taking two steps gives same estimates, stick with the first step.
        if np.isnan(self.var).any():
            self.var = self.oldvar[1]

        self.xi = np.exp(digamma(self.var)-digamma(insum(self.var,[1])))
        self.require()

