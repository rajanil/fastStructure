
import numpy as np
cimport numpy as np
cimport admixprop as ap
from cpython cimport bool
from scipy.special import digamma, gammaln, polygamma
import scipy.optimize as opt
import utils
from functools import reduce

ctypedef np.uint8_t uint8_t

cdef extern from "allelefreq.h":
    void P_update_simple( uint8_t* G, double* zetabeta, double* zetagamma, double* xi, double* beta, double* gamma, double* var_beta, double* var_gamma, long N, long L, long K )
    void P_update_logistic( double* Dvarbeta, double* Dvargamma, double* mu, double* Lambda, double* var_beta, double* var_gamma, double mintol, long L, long K)

cdef class AlleleFreq:

    def __cinit__(self, long L, long K, str prior):

        """
        Sets initial parameter values for the variational distributions over
        allele frequencies. 
        """

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
        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.oldvar_beta = []
        self.oldvar_gamma = []
        self.require()

    cdef copy(self):

        """
        Creates a new instance of the class with explicit copies of relevant
        variables.
        """

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

        """
        Enforces variables of type `numpy.ndarray` to be in C-contiguous order.
        """

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

    cdef _update_simple(self, np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi):

        """
        Update parameters of variational distributions over
        allele frequencies, given genotype data and estimates
        of parameters of variational distributions over admixture 
        proportions. This update method is called when the
        ``simple prior'' over allele frequencies is chosen.

        Arguments
            G : numpy array of genotypes

            psi : instance of `AdmixProp`

        """

        self.var_beta = np.zeros((self.L,self.K),dtype=np.float64)
        self.var_gamma = np.zeros((self.L,self.K),dtype=np.float64)
        self.require()

        P_update_simple(<np.uint8_t*> G.data, <double*> self.zetabeta.data, <double*> self.zetagamma.data, <double*> psi.xi.data, <double*> self.beta.data, <double*> self.gamma.data, <double*> self.var_beta.data, <double*> self.var_gamma.data, psi.N, self.L, self.K)

        # if the update fails for some reason, stick with the previous set of values
        if np.isnan(self.var_beta).any():
            self.var_beta = self.oldvar_beta[-1]

        if np.isnan(self.var_gamma).any():
            self.var_gamma = self.oldvar_gamma[-1]

        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.require()

    cdef _update_logistic(self, np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi):

        """
        Update parameters of variational distributions over
        allele frequencies, given genotype data and estimates
        of parameters of variational distributions over admixture
        proportions. This update method is called when the
        ``logistic prior'' over allele frequencies is chosen.

        Arguments
            G : numpy array of genotypes

            psi : instance of `AdmixProp`

        """

        cdef np.ndarray Dvarbeta, Dvargamma, beta, var_beta, var_gamma, vars_at_boundary, varbeta, vargamma, bad_beta, bad_gamma
        cdef list bad_conditions
        beta = np.zeros((self.L,self.K),dtype=np.float64)

        # compute data-dependent terms in the update equations
        Dvarbeta = self.var_beta.copy()
        Dvarbeta = np.require(Dvarbeta, dtype=np.float64, requirements='C')
        Dvargamma = self.var_gamma.copy()
        Dvargamma = np.require(Dvargamma, dtype=np.float64, requirements='C')

        P_update_simple(<np.uint8_t*> G.data, <double*> self.zetabeta.data, <double*> self.zetagamma.data, <double*> psi.xi.data, <double*> beta.data, <double*> beta.data, <double*> Dvarbeta.data, <double*> Dvargamma.data, psi.N, self.L, self.K)

        # use an iterative fixed-point solver to update parameter estimates.
        var_beta, var_gamma = self._unconstrained_solver(Dvarbeta, Dvargamma)

        # if a variable violates positivity constraint, 
        # set it to an estimate from the previous update
        bad_beta = reduce(utils.OR,[(var_beta<=0),np.isnan(var_beta)])
        bad_gamma = reduce(utils.OR,[(var_gamma<=0),np.isnan(var_gamma)])
        var_beta[bad_beta] = self.var_beta[bad_beta]
        var_gamma[bad_gamma] = self.var_gamma[bad_gamma]

        self.var_beta = var_beta
        self.var_gamma = var_gamma
        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.require()

    cdef _unconstrained_solver(self, np.ndarray[np.float64_t, ndim=2] Dvarbeta, np.ndarray[np.float64_t, ndim=2] Dvargamma):

        """
        Iterative fixed-point solver to update estimates of variational parameters
        for allele frequencies, when the logistic prior is chosen. The update equations
        in this case, have the same form as those for the simple prior.

        Arguments:

            Dvarbeta : numpy.ndarray
                data-dependent terms relevant for the update of `beta` parameters 

            Dvargamma : numpy.ndarray
                data-dependent terms relevant for the update of `gamma` parameters

        Note: 
            positivity constraints on variables are not explicitly
            enforced in this iterative scheme.
        """

        cdef np.ndarray var_beta, var_gamma
        var_beta = self.var_beta.copy()
        var_beta = np.require(var_beta, dtype=np.float64, requirements='C')
        var_gamma = self.var_gamma.copy()
        var_gamma = np.require(var_gamma, dtype=np.float64, requirements='C')

        P_update_logistic(<double*> Dvarbeta.data, <double*> Dvargamma.data, <double*> self.mu.data, <double*> self.Lambda.data, <double*> var_beta.data, <double*> var_gamma.data, self.mintol, self.L, self.K)

        return var_beta, var_gamma

    cdef update(self, np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi):

        """
        Calls the relevant update method depending on the choice of prior.

        Arguments
            G : numpy array of genotypes

            pi : instance of `AlleleFreq`

        """

        if self.prior=='simple':
            self._update_simple(G, psi)
        elif self.prior=='logistic':
            self._update_logistic(G, psi)

    cdef square_update(self, np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi):

        """
        Accelerated update of variational parameters of
        allele frequencies.

        Arguments
            G : numpy array of genotypes

            pi : instance of `AlleleFreq`

        """

        cdef long step
        cdef bool a_ok
        cdef np.ndarray R_beta, R_gamma, V_beta, V_gamma
        cdef double a

        self.oldvar_beta = [self.var_beta.copy()]
        self.oldvar_gamma = [self.var_gamma.copy()]

        # take two update steps
        for step from 0 <= step < 2:
            self.update(G, psi)
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

        # given two update steps, compute an optimal step that achieves
        # a better marginal likelihood than the best of the two steps.
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

        # if this accelerated step fails for some reason, stick with the first non-accelerated step.
        if np.isnan(self.var_beta).any() or np.isnan(self.var_gamma).any():
            self.var_beta = self.oldvar_beta[1]
            self.var_gamma = self.oldvar_gamma[1]

        self.zetabeta = np.exp(digamma(self.var_beta) - digamma(self.var_beta+self.var_gamma))
        self.zetagamma = np.exp(digamma(self.var_gamma) - digamma(self.var_beta+self.var_gamma))
        self.require()

    cdef update_hyperparam(self, bool nolambda):

        """
        Update parameters of the logistic prior over allele frequencies.

        Arguments:

            nolambda : bool
                if True, the Lambda hyperparameter is NOT updated.

        """

        cdef np.ndarray dat, C
        if self.prior=='logistic':
            dat = digamma(self.var_beta)-digamma(self.var_gamma)
            self.mu = utils.insum(self.Lambda*dat,[1]) / self.Lambda.sum()
            diff = dat-self.mu

            if not nolambda:
                C = 1./(self.L) * (utils.outsum(diff**2) + utils.outsum(polygamma(1,self.var_beta)+polygamma(1,self.var_gamma))).ravel()
                self.Lambda = 1./C
