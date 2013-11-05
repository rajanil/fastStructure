import numpy as np
cimport numpy as np
cimport allelefreq as af
cimport admixprop as ap
from scipy.special import digamma, gammaln, polygamma
import utils

ctypedef np.uint8_t uint8_t

cdef extern from "marglikehood.h":
    double marglikehood( uint8_t* G, double* zetabeta, double* zetagamma, double* xi, long N, long L, long K )

cdef double marginal_likelihood(np.ndarray[np.uint8_t, ndim=2] G, ap.AdmixProp psi, af.AlleleFreq pi):

    cdef double E1, E2, E3, Etotal

    E1 = marglikehood(<np.uint8_t*> G.data, <double*> pi.zetabeta.data, <double*> pi.zetagamma.data, <double*> psi.xi.data, psi.N, pi.L, pi.K)

    E2 = (utils.insum(gammaln(psi.var) - gammaln(psi.alpha) - (psi.var-psi.alpha)*np.nan_to_num(np.log(psi.xi)),[1]) \
        - gammaln(utils.insum(psi.var,[1])) + gammaln(utils.insum(psi.alpha,[1]))).sum()

    if pi.prior=='simple':
        E3 = (gammaln(pi.var_beta) - gammaln(pi.beta) - (pi.var_beta-pi.beta)*np.nan_to_num(np.log(pi.zetabeta)) \
            + gammaln(pi.var_gamma) - gammaln(pi.gamma) - (pi.var_gamma-pi.gamma)*np.nan_to_num(np.log(pi.zetagamma)) \
            - gammaln(pi.var_beta+pi.var_gamma) + gammaln(pi.beta+pi.gamma)).sum()
    elif pi.prior=='logistic':
        diff = digamma(pi.var_beta)-digamma(pi.var_gamma)-pi.mu
        E3 = 0.5*pi.L*np.log(pi.Lambda).sum() - 0.5*(pi.Lambda*utils.outsum(diff**2)).sum() \
            - 0.5*(pi.Lambda*utils.outsum(polygamma(1,pi.var_beta)+polygamma(1,pi.var_gamma))).sum() \
            - np.sum(utils.nplog(pi.zetabeta)+utils.nplog(pi.zetagamma)) \
            + ((pi.var_beta>0)*gammaln(pi.var_beta) - (pi.var_beta-1)*utils.nplog(pi.zetabeta) \
            + (pi.var_gamma>0)*gammaln(pi.var_gamma) - (pi.var_gamma-1)*utils.nplog(pi.zetagamma) \
            - gammaln(pi.var_beta+pi.var_gamma)).sum()

    Etotal = (E1 + E2 + E3)/float(psi.N*pi.L)

    return Etotal
