
import numpy as np
cimport numpy as np
from cpython cimport bool
import vars.utils as utils
cimport vars.allelefreq as af 
import vars.allelefreq as af
cimport vars.admixprop as ap 
import vars.admixprop as ap
import vars.marglikehood as mlhood
import time

def infer_variational_parameters(np.ndarray[np.uint8_t, ndim=2] G, int K, str outfile, double mintol, str prior, int cv):

    cdef int N, L, batch_size, restart, iter
    cdef double Estart, E, Enew, reltol, diff, totaltime, itertime
    cdef np.ndarray g
    cdef list indices, old, Psis, Es, Times
    cdef ap.AdmixProp psi, psistart
    cdef af.AlleleFreq pi, pistart, piG

    totaltime = time.time()
    N = G.shape[0]
    L = G.shape[1]

    # minimum number of SNPs that can resolve populations
    # with an Fst of 0.0001, given sample size `N`
    batch_size = min([L,int(1000000/N)])

    handle = open('%s.%d.log'%(outfile,K),'w')
    handle.close()

    itertime = time.time()

    # First, select initial values for the variational parameters
    # from 5 random initializations of the algorithm.
    Estart = -np.inf
    for restart in xrange(5):

        # if dataset is too large, initialize variational parameters 
        # using a random subset of the data (similar to a warm start).
        if batch_size<L:

            # choose random subset of the SNPs
            indices = list(utils.random_combination(xrange(L), batch_size))
            g = G[:,indices]
            g = np.require(g, dtype=np.uint8, requirements='C')

            # iterate the variational algorithm to `weak' convergence
            # to initialize parameters for admixture proportions
            psi = ap.AdmixProp(N, K)
            pi = af.AlleleFreq(batch_size, K, prior)
            # variational admixture proportion update
            psi.update(g, pi)
            # variational allele frequency update
            pi.update(g, psi)
            E = mlhood.marginal_likelihood(g, psi, pi)
            reltol = np.inf
            iter = 0
            while np.abs(reltol)>mintol*10:
                # accelerated variational admixture proportion update
                psi.square_update(g, pi)
                # accelerated variational allele frequency update
                pi.square_update(g, psi)
                iter += 1
                if iter%10==0:
                    Enew = mlhood.marginal_likelihood(g, psi, pi)
                    reltol = Enew - E
                    E = Enew
                    # update allele frequency hyperparameters
                    pi.update_hyperparam(False)
            piG = pi.copy()

            # after initializing admixture proportions, initialize the allele
            # frequency parameters for all SNPs, keeping admixture proportions fixed.
            # if the logistic prior is chosen, hyperparameter `Lambda`
            # is also kept fixed in this step.
            pi = af.AlleleFreq(L, K, prior)
            if pi.prior=='logistic':
                pi.Lambda = piG.Lambda.copy()
            old = [pi.var_beta.copy(), pi.var_gamma.copy()]
            diff = np.inf
            while diff>1e-1:
                # accelerated variational allele frequency update
                pi.square_update(G, psi)
                diff = np.mean(np.abs(pi.var_beta-old[0])+np.abs(pi.var_gamma-old[1]))
                old = [pi.var_beta.copy(), pi.var_gamma.copy()]
                # update allele frequency hyperparameters
                pi.update_hyperparam(True)

        else:

            # simple initializing of variational parameters (cold start)
            psi = ap.AdmixProp(N,K)
            pi = af.AlleleFreq(L, K, prior)

        # compute marginal likelihood for this initialization
        psi.update(G, pi)
        pi.update(G, psi)
        E = mlhood.marginal_likelihood(G, psi, pi)
        handle = open('%s.%d.log'%(outfile,K),'a')
        handle.write("Marginal likelihood with initialization (%d) = %.10f\n"%(restart+1,E))
        handle.close()
        # select current initialization if it has a higher marginal
        # likelihood than all previous initializations
        if E>Estart:
            pistart = pi.copy()
            psistart = psi.copy()
            Estart = E

    itertime = time.time()-itertime
    E = Estart
    pi = pistart.copy()
    psi = psistart.copy()
    iter = 0

    handle = open('%s.%d.log'%(outfile,K),'a')
    to_write = ["Iteration", "Marginal_Likelihood", "delta_Marginal_Likelihood", "Iteration_Time (secs)"]
    handle.write(' '.join(to_write)+"\n")
    to_write = ['%d'%iter, '%.10f'%E, '--', '%.3f'%itertime]
    handle.write(' '.join(to_write)+"\n")
    handle.close()

    itertime = time.time()
    reltol = np.inf
    while np.abs(reltol)>mintol:
        # accelerated variational admixture proportion update
        psi.square_update(G, pi)
        psi.update(G, pi)

        # accelearted variational allele frequency update
        pi.square_update(G, psi)
        pi.update(G, psi)

        # Compute marginal likelihood once every 10 iterations
        if (iter+1)%10==0:

            E_new = mlhood.marginal_likelihood(G, psi, pi)
            reltol = E_new-E
            E = E_new
            itertime = time.time()-itertime

            handle = open('%s.%d.log'%(outfile,K),'a')
            to_write = ['%d'%(iter+1), '%.10f'%E, '%.10f'%reltol, '%.3f'%itertime]
            handle.write(' '.join(to_write)+"\n")
            handle.close()

            itertime = time.time()
            pi.update_hyperparam(False)

        iter += 1

    # posterior mean of allele frequencies and admixture proportions
    P = pi.var_beta/(pi.var_beta+pi.var_gamma)
    Q = psi.var/utils.insum(psi.var,[1])

    totaltime = time.time()-totaltime
    handle = open('%s.%d.log'%(outfile,K),'a')
    handle.write("Marginal Likelihood = %.10f\n"%E)
    handle.write("Total time = %.4f seconds\n"%totaltime)
    handle.write("Total iterations = %d \n"%iter)
    handle.close()

    # Computing cross-validation error
    if cv:
        meandeviance = CV(G, psi, pi, cv, mintol)
        handle = open('%s.%d.log'%(outfile,K),'a')
        handle.write("CV error = %.7f, %.7f \n"%(np.mean(meandeviance), np.std(meandeviance, ddof=1)))
        handle.close()

    other = dict([('varQ',psi.var), ('varPb',pi.var_beta), ('varPg',pi.var_gamma)])

    return Q, P, other

cdef double expected_genotype(ap.AdmixProp psi, af.AlleleFreq pi, int n, int l):

    """
    compute the expected genotype of missing data, given observed genotypes, 
    after integrating out latent variables and model parameters using
    their variational distributions.

    Arguments:

        psi : instance of `AdmixProp`
            variational distribution of admixture proprotions
        
        pi : instance of `AlleleFreq`
            variational distribution of allele frequencies

        n : int
            sample index

        l : int
            SNP index

    Returns:

        float

    """

    cdef np.ndarray g = np.zeros((3,),dtype=float)
    cdef np.ndarray Q, Qi, Pb, Pib, Pg, Pig, P
    cdef double nu

    # selecting relevant parameters
    Q = psi.xi[n:n+1]
    Qi = Q/Q.sum()
    Pb = pi.var_beta[l:l+1]
    Pg = pi.var_gamma[l:l+1]
    Pib = Pb/(Pb+Pg)
    Pig = Pg/(Pb+Pg)

    # probability that genotype is zero
    P = Pig.T*Pig
    P[range(pi.K),range(pi.K)] = (Pig*(Pg+1)/(Pb+Pg+1)).ravel()
    P = np.dot(Qi,np.dot(P,Qi.T))
    g[0] = np.sum(P)

    # probability that genotype is one
    P = Pib.T*Pig
    P[range(pi.K),range(pi.K)] = (Pib*Pg/(Pb+Pg+1)).ravel()
    P = np.dot(Qi,np.dot(P,Qi.T))
    g[1] = 2*np.sum(P)

    # probability that genotype is two
    P = Pib.T*Pib
    P[range(pi.K),range(pi.K)] = (Pib*(Pb+1)/(Pb+Pg+1)).ravel()
    P = np.dot(Qi,np.dot(P,Qi.T))
    g[2] = np.sum(P)

    nu = np.sum(g*np.arange(3))

    return nu

cdef np.ndarray CV(np.ndarray[np.uint8_t, ndim=2] Gtrue, ap.AdmixProp psi, af.AlleleFreq pi, int cv, double mintol):

    """
    compute the cross-validation error for a dataset, by computing the
    model deviance on held-out subsets of the data. 

    Arguments:

        Gtrue : numpy.ndarray
            array of genotypes
    
        psi : instance of `admixprop.AdmixProp`
            seed for variational admixture proportion parameters

        pi : instance of `allelefreq.AlleleFreq`
            seed for variational allele frequency parameters

        cv : int
            number of folds of cross-validation

        mintol : double
            convergence criterion for the variational algorithm

    Returns:

        numpy.ndarray

    Note:
        Given the variational parameters estimated using the entire
        dataset, a random subset of the genotype entries are held-out
        and the variational algorithm is run on the remaining genotypes,
        with the estimated parameters as the initial seed. This
        allows for fast parameter estimation. The prediction error
        for the held-out dataset is computed from the binomial deviance
        of the held-out genotype entries, given the re-estimated 
        variational parameters.

    """

    cdef bool wellmasked = False
    cdef list nonindices, masks, newmasks, deviances
    cdef tuple mask
    cdef int i, j, n, l, m, iter, N, L
    cdef double E, E_new, reltol, deviance
    cdef np.ndarray G, Gmask, pg, meandeviance
    cdef ap.AdmixProp psimask
    cdef af.AlleleFreq pimask

    N = Gtrue.shape[0]
    L = Gtrue.shape[1]
    while not wellmasked:
        # partition the observed genotype entries
        nonindices = [(i,j) for i,j in zip((Gtrue==3).nonzero()[0], (Gtrue==3).nonzero()[1])]
        masks = utils.random_partition([Gtrue.shape[0], Gtrue.shape[1]], nonindices, cv)
        wellmasked = True

        # test to ensure that for all partitions, the loci are all variant
        newmasks = []
        for mask in masks:
            G = Gtrue.copy()
            Gmask = -1*np.ones((N,L), dtype='int8')
            Gmask[mask[0],mask[1]] = G[mask[0],mask[1]]
            G[mask[0],mask[1]] = 3
            if not (((G==1)+(G==2)).sum(0)==0).any():
                newmasks.append(mask)

        if not len(newmasks)>=cv:
            wellmasked = False
            print "Failed"

    masks = newmasks[:cv]
    meandeviance = np.zeros((cv,), dtype=float)
    for m,mask in enumerate(masks):
        deviances = []
        # construct a masked genotype matrix
        G = Gtrue.copy()
        Gmask = -1*np.ones((N,L), dtype='int8')
        Gmask[mask[0],mask[1]] = G[mask[0],mask[1]]
        G[mask[0],mask[1]] = 3

        # restimate parameters from masked genotype matrix
        psimask = psi.copy()
        pimask = pi.copy()
        E = mlhood.marginal_likelihood(G, psi, pi)
        reltol = np.inf
        iter = 0

        while np.abs(reltol)>mintol:
            # VBE step
            psimask.square_update(G, pimask)
            psimask.update(G, pimask)

            # VBM step
            pimask.square_update(G, psimask)
            pimask.update(G, psimask)

            # optimize hyperparameters once every 10 iterations
            if (iter+1)%10==0:

                E_new = mlhood.marginal_likelihood(G, psimask, pimask)
                reltol = E_new-E
                E = E_new
                pimask.update_hyperparam(False)

            iter += 1

        for n,gmask in enumerate(Gmask):
            for l in (gmask!=-1).nonzero()[0]:
                nu = expected_genotype(psimask, pimask, n, l)
                deviance = gmask[l]*utils.nplog(gmask[l]/nu) + (2-gmask[l])*utils.nplog((2-gmask[l])/(2-nu))
                deviances.append(deviance)
        meandeviance[m] = np.mean(deviances)

    return meandeviance
