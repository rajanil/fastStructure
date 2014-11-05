import numpy as np
cimport numpy as np
cimport allelefreq as af
from scipy.special import digamma
from cpython cimport bool
import utils

ctypedef np.uint8_t uint8_t

cdef extern from "admixprop.h":
    void Q_update( uint8_t* G, double* zetabeta, double* zetagamma, double* xi, double* new_var, long N, long L, long K )

cdef class AdmixProp:

    def __cinit__(self, long N, long K):

        """
        Sets initial parameter values for the variational distributions over
        admixture proportions. The prior over admixture proportions is set 
        to be a symmetric Dirichlet distribution with parameter 1/K.
        """

        self.N = N
        self.K = K

        # Initializing hyperparameters
        self.alpha = 1./K*np.ones((1,K))

        # Initializing variational parameters for admixture proportions
        self.var = np.ones((N,K)) + 0.1*np.random.rand(N,K)
        self.xi = np.exp(digamma(self.var)-digamma(utils.insum(self.var,[1])))
        self.oldvar = []

    cdef copy(self):

        """
        Creates a new instance of the class with explicit copies of relevant
        variables.
        """

        cdef AdmixProp newinstance
        newinstance = AdmixProp(self.N, self.K)
        newinstance.var = self.var.copy()
        newinstance.xi = self.xi.copy()
        newinstance.oldvar = []

        return newinstance

    cdef require(self):

        """
        Enforces variables of type `numpy.ndarray` to be in C-contiguous order.
        """

        self.var = np.require(self.var, dtype=np.float64, requirements='C')
        self.xi = np.require(self.xi, dtype=np.float64, requirements='C')

    cdef update(self, np.ndarray[np.uint8_t, ndim=2] G, af.AlleleFreq pi):

        """
        Update parameters of variational distributions over 
        admixture proportions, given genotype data and estimates
        of parameters of variational distributions over allele
        frequencies.

        Arguments
            G : numpy array of genotypes

            pi : instance of `AlleleFreq`

        """

        self.var = np.zeros((self.N,self.K), dtype=np.float64)
        Q_update(<np.uint8_t*>G.data, <double*> pi.zetabeta.data, <double*> pi.zetagamma.data, <double*> self.xi.data, <double*> self.var.data, self.N, pi. L, self.K)

        # if the update fails for some reason, stick with the previous set of values
        if np.isnan(self.var).any():
            self.var = self.oldvar[-1]
        else:
            self.var = self.alpha + self.var
            self.xi = np.exp(digamma(self.var)-digamma(utils.insum(self.var,[1])))
        self.require()

    cdef square_update(self, np.ndarray[np.uint8_t, ndim=2] G, af.AlleleFreq pi):

        """
        Accelerated update of variational parameters of 
        admixture proportions.

        Arguments
            G : numpy array of genotypes

            pi : instance of `AlleleFreq`

        """

        cdef long step
        cdef bool a_ok
        cdef np.ndarray R, V
        cdef double a

        self.oldvar = [self.var.copy()]
        # take two update steps
        for step from 0 <= step < 2:
            self.update(G, pi)
            self.oldvar.append(self.var.copy())

        R = self.oldvar[1] - self.oldvar[0]
        V = self.oldvar[2] - self.oldvar[1] - R
        a = -1.*np.sqrt((R*R).sum()/(V*V).sum())

        if a>-1:
            a = -1.

        # given two update steps, compute an optimal step that achieves
        # a better marginal likelihood than the best of the two steps.
        a_ok = False
        while not a_ok:
            self.var = (1+a)**2*self.oldvar[0] - 2*a*(1+a)*self.oldvar[1] + a**2*self.oldvar[2]
            if (self.var<=0).any():
                a = (a-1)/2.
                if np.abs(a+1)<1e-4:
                    a = -1.
            else:
                a_ok = True

        # if this accelerated step fails for some reason, stick with the first non-accelerated step
        if np.isnan(self.var).any():
            self.var = self.oldvar[1]

        self.xi = np.exp(digamma(self.var)-digamma(utils.insum(self.var,[1])))
        self.require()

