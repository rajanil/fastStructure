
import numpy as np
cimport numpy as np
from cpython cimport bool
import ctypes
from scipy.special import digamma, gammaln, polygamma
import scipy.optimize as opt
import utils

ctypedef np.uint8_t uint8_t

cdef extern from "allelefreq.h":
    void P_update_simple( uint8_t* G, double* zetabeta, double* zetagamma, double* xi, double* beta, double* gamma, double* var_beta, double* var_gamma, long N, long L, long K )
    void P_update_logistic( double* Dvarbeta, double* Dvargamma, double* mu, double* Lambda, double* var_beta, double* var_gamma, double mintol, long L, long K)

cdef class AlleleFreq:
    """Parameter of variational distribution over genotype frequencies.

    Arguments
        L : int

        K : int

    """

    cdef long L,K
    cdef double mintol
    cdef np.ndarray beta, gamma, var_beta, var_gamma, zetabeta, zetagamma, piA, F, mu, Lambda
    cdef list oldvar_beta, oldvar_gamma
    cdef str prior 

    def __cinit__(self, long L, long K, np.ndarray[np.uint8_t, ndim=2] G=np.empty((1,1),dtype=np.uint8), str prior='simple'):

        self.L = L
        self.K = K
        self.prior = prior
        if self.prior=='simple':
            self.beta = np.ones((self.L,self.K))
            self.gamma = np.ones((self.L,self.K))
        elif self.prior=='logistic':
            self.mu = np.zeros((self.L,1))
            self.Lambda = np.ones((self.K,))
            self.mintol = 1e-1

        self.var_beta = np.ones((L,K)) + 0.1*np.random.rand(L,K)
        self.var_gamma = 10*np.ones((L,K)) + 0.1*np.random.rand(L,K)
#        self.var_beta = np.ones((L,K))
#        self.var_gamma = 10*np.ones((L,K))
        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.oldvar_beta = []
        self.oldvar_gamma = []
        self.require()

    cdef copy(self):

        cdef AlleleFreq newinstance
        newinstance = AlleleFreq(self.L, self.K, prior=self.prior)
        newinstance.var_beta = self.var_beta.copy()
        newinstance.zetabeta = self.zetabeta.copy()
        newinstance.var_gamma = self.var_gamma.copy()
        newinstance.zetagamma = self.zetagamma.copy()

        if self.prior=='logistic':
            newinstance.mu = self.mu
            newinstance.Lambda = self.Lambda

        newinstance.require()
        return newinstance

    cdef require(self):

        self.var_beta = np.require(self.var_beta, dtype=np.float64, requirements='C')
        self.var_gamma = np.require(self.var_gamma, dtype=np.float64, requirements='C')
        self.zetabeta = np.require(self.zetabeta, dtype=np.float64, requirements='C')
        self.zetagamma = np.require(self.zetagamma, dtype=np.float64, requirements='C')
        if self.prior=='simple':
            self.beta = np.require(self.beta, dtype=np.float64, requirements='C')
            self.gamma = np.require(self.gamma, dtype=np.float64, requirements='C')
        elif self.prior=='logistic':
            self.mu = np.require(self.mu, dtype=np.float64, requirements='C')
            self.Lambda = np.require(self.Lambda, dtype=np.float64, requirements='C')

    cdef _update_simple(self, np.ndarray[np.uint8_t, ndim=2] G, Psi psi):
        """Update parameters of distribution over allele frequencies in the VBM step.

        Arguments
            G : array

        """

        self.var_beta = np.zeros((self.L,self.K),dtype=np.float64)
        self.var_gamma = np.zeros((self.L,self.K),dtype=np.float64)
        self.require()

        P_update(<np.uint8_t*> G.data, <double*> self.zetabeta.data, <double*> self.zetagamma.data, <double*> psi.xi.data, <double*> self.beta.data, <double*> self.gamma.data, <double*> self.var_beta.data, <double*> self.var_gamma.data, psi.N, self.L, self.K)

        if np.isnan(self.var_beta).any():
            self.var_beta = self.oldvar_beta[-1]

        if np.isnan(self.var_gamma).any():
            self.var_gamma = self.oldvar_gamma[-1]

        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.require()

    cdef _update_logistic(self, np.ndarray[np.uint8_t, ndim=2] G, Psi psi):

        # compute contribution from data
        cdef np.ndarray Dvarbeta, Dvargamma, beta, var_beta, var_gamma, vars_at_boundary, varbeta, vargamma, bad_beta, bad_gamma
        cdef list bad_conditions
        beta = np.zeros((self.L,self.K),dtype=np.float64)

        # compute contribution from data
        Dvarbeta = self.var_beta.copy()
        Dvarbeta = np.require(Dvarbeta, dtype=np.float64, requirements='C')
        Dvargamma = self.var_gamma.copy()
        Dvargamma = np.require(Dvargamma, dtype=np.float64, requirements='C')

        P_update(<np.uint8_t*> G.data, <double*> self.zetabeta.data, <double*> self.zetagamma.data, <double*> psi.xi.data, <double*> beta.data, <double*> beta.data, <double*> Dvarbeta.data, <double*> Dvargamma.data, psi.N, self.L, self.K)

        # using iterative solver
        var_beta, var_gamma = self._unconstrained_solver(Dvarbeta, Dvargamma)

        # using constrained optimiation solver
        bad_conditions = [(var_beta<=0),(var_gamma<=0),np.isnan(var_beta),np.isnan(var_gamma)]
        vars_at_boundary = np.any(reduce(OR,bad_conditions),1)

        if vars_at_boundary.sum():
            print "%d vars at boundary"%vars_at_boundary.sum()
#            var_beta[vars_at_boundary] = self.var_beta[vars_at_boundary]
#            var_gamma[vars_at_boundary] = self.var_gamma[vars_at_boundary]
#            varbeta, vargamma = self._constrained_solver(G, Dvarbeta, Dvargamma, vars_at_boundary)
#            var_beta[vars_at_boundary,:] = varbeta.copy()
#            var_gamma[vars_at_boundary,:] = vargamma.copy()

        bad_beta = reduce(OR,[(var_beta<=0),np.isnan(var_beta)])
        bad_gamma = reduce(OR,[(var_gamma<=0),np.isnan(var_gamma)])
        var_beta[bad_beta] = self.var_beta[bad_beta]
        var_gamma[bad_gamma] = self.var_gamma[bad_gamma]

        self.var_beta = var_beta
        self.var_gamma = var_gamma
        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.require()

    cdef _unconstrained_solver(self, np.ndarray[np.float64_t, ndim=2] Dvarbeta, np.ndarray[np.float64_t, ndim=2] Dvargamma):

        cdef np.ndarray var_beta, var_gamma
        var_beta = self.var_beta.copy()
        var_beta = np.require(var_beta, dtype=np.float64, requirements='C')
        var_gamma = self.var_gamma.copy()
        var_gamma = np.require(var_gamma, dtype=np.float64, requirements='C')

        P_update_logistic(<double*> Dvarbeta.data, <double*> Dvargamma.data, <double*> self.mu.data, <double*> self.Lambda.data, <double*> var_beta.data, <double*> var_gamma.data, self.mintol, self.L, self.K)

        return var_beta, var_gamma

    cdef _constrained_solver(self, np.ndarray[np.uint8_t, ndim=2] G, np.ndarray[np.float64_t, ndim=2] Dvarbeta, np.ndarray[np.float64_t, ndim=2] Dvargamma, constraint):

        cdef long v, L, numvars
        cdef np.ndarray varbeta, vargamma, zetabeta, zetagamma, dvarbeta, dvargamma, mu, Lambda
        cdef list bounds
        cdef np.ndarray vbeta, vgamma, lgbeta, lggamma, lgbetagamma, dbeta, dgamma, dbetagamma, pbeta, pgamma, pbetagamma, ppbeta, ppgamma, diff, Dfbeta, Dfgamma, Df, func, xo
        cdef double f

        L = constraint.sum()
        varbeta = self.var_beta[constraint,:]
        vargamma = self.var_gamma[constraint,:]
        zetabeta = self.zetabeta[constraint,:]
        zetagamma = self.zetagamma[constraint,:]
        dvarbeta = Dvarbeta[constraint,:]
        dvargamma = Dvargamma[constraint,:]
        mu = self.mu[constraint]

        numvars = 2*L*self.K
        bounds = [(0,np.inf) for v in xrange(numvars)]

        return varbeta, vargamma

    cdef update(self, np.ndarray[np.uint8_t, ndim=2] G, Psi psi):

        if self.prior=='simple':
            self._update_simple(G, psi)
        elif self.prior=='logistic':
            self._update_logistic(G, psi)

    cdef square_update(self, np.ndarray[np.uint8_t, ndim=2] G, Psi psi):

        cdef long step
        cdef bool a_ok
        cdef np.ndarray R_beta, R_gamma, V_beta, V_gamma
        cdef double a

        self.oldvar_beta = [self.var_beta.copy()]
        self.oldvar_gamma = [self.var_gamma.copy()]

        for step from 0 <= step < 2:
            self.update_VBM(G, psi)
            self.oldvar_beta.append(self.var_beta.copy())
            self.oldvar_gamma.append(self.var_gamma.copy())

        R_beta = self.oldvar_beta[1] - self.oldvar_beta[0]
        R_gamma = self.oldvar_gamma[1] - self.oldvar_gamma[0]
        V_beta = self.oldvar_beta[2] - self.oldvar_beta[1] - R_beta
        V_gamma = self.oldvar_gamma[2] - self.oldvar_gamma[1] - R_gamma

        a = -1.*np.sqrt(((R_beta*R_beta).sum()+(R_gamma*R_gamma).sum())
                / ((V_beta*V_beta).sum()+(V_gamma*V_gamma).sum()))

        if a>-1:
            a = -1.

        a_ok = False
        while not a_ok:
            self.var_beta = (1+a)**2*self.oldvar_beta[0] - 2*a*(1+a)*self.oldvar_beta[1] + a**2*self.oldvar_beta[2]
            self.var_gamma = (1+a)**2*self.oldvar_gamma[0] - 2*a*(1+a)*self.oldvar_gamma[1] + a**2*self.oldvar_gamma[2]
            if (self.var_beta<=0).any() or (self.var_gamma<=0).any():
                a = (a-1)/2.
                if np.abs(a+1)<1e-4:
                    a = -1.
            else:
                a_ok = True

        if np.isnan(self.var_beta).any() or np.isnan(self.var_gamma).any():
            self.var_beta = self.oldvar_beta[1]
            self.var_gamma = self.oldvar_gamma[1]

        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.require()

    cdef update_hyperparam(self, bool nolambda):

        cdef np.ndarray dat, C
        if self.prior=='logistic':
            dat = digamma(self.var_beta)-digamma(self.var_gamma)
            self.mu = insum(self.Lambda*dat,[1]) / self.Lambda.sum()
            diff = dat-self.mu

            if not nolambda:
                C = 1./(self.L) * (outsum(diff**2) + outsum(polygamma(1,self.var_beta)+polygamma(1,self.var_gamma))).ravel()
                self.Lambda = 1./C
