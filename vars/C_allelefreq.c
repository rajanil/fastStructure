
#include "allelefreq.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_errno.h>

void P_update_simple(const uint8_t* G, const double* zetabeta, const double* zetagamma, const double* xi, const double* beta, const double* gamma, double* var_beta, double* var_gamma, long N, long L, long K)
{
    uint8_t genotype;
    long idx, n, l, k;
    double theta_beta_sum, theta_gamma_sum;
    double *var_beta_tmp, *var_gamma_tmp;

    var_beta_tmp = (double*) malloc(K * sizeof(double));
    var_gamma_tmp = (double*) malloc(K * sizeof(double));

    // loop over loci
    for (l=0; l<L; l++) {

        for (k=0; k<K; k++) {

            var_beta_tmp[k] = 0.0;
            var_gamma_tmp[k] = 0.0;
        }

        // loop over samples
        for (n=0; n<N; n++) {

            genotype = G[n*L+l];

            // missing data do not contribute
            if (genotype!=3) {

	        // Add across the components to compute the normalizing constant
                // for the indicator variables.
                // These are already exponentiated, so they are terms from exp(log likelihood).
                // xi: exp(E(log(probability of population k for individual n))
                // zetabeta: exp(E(log(probability of allele l for population k))
                // compute xi*zeta_{beta,gamma}
                theta_beta_sum = 0.0;
                theta_gamma_sum = 0.0;
                for (k=0; k<K; k++) {
		    // In my notation:
		    // theta_beta is the indicator that allele A is in population k.
		    // theta_gamma is the indicator that allele A is in population k.
                    //theta_beta_sum += xi[n * K + k] * zetabeta[l * K + k];
                    //theta_gamma_sum += xi[n * K + k] * zetagamma[l * K + k];
   		                        theta_beta_sum += (double) genotype * xi[n * K + k] * zetabeta[l * K + k];
                    theta_gamma_sum += (double) (2 - genotype) * xi[n * K + k] * zetagamma[l * K + k];
                }
                theta_beta_sum = theta_beta_sum == 0.0 ? 1.0: theta_beta_sum;
                theta_gamma_sum = theta_gamma_sum == 0.0 ? 1.0: theta_gamma_sum;

                // increment var_{beta,gamma}_tmp
                for (k=0; k<K; k++) {
		  // genotype is either 0, 1, or 2.
                  // If it is 2, both alleles count towards beta.  
                  // If it is 0, both alleles count towards gamma.
                  // If it is 1, each allele gets one.
                  // Note that this is all multiplied by zetabeta and zetagamma below.
                  //
                  // TOOD: I believe this is a bug.
                  // Why do the denominators sum over all the terms, but the numerators
                  // only by the genotypes?  These indicators do not sum to one?
                    var_beta_tmp[k] += (double) genotype * xi[n * K + k] / theta_beta_sum;
                    var_gamma_tmp[k] += (double) (2 - genotype) * xi[n * K + k] / theta_gamma_sum;
                }
            }
        }

        // compute var_{beta,gamma}
        for (k=0; k<K; k++) {
	  // The variables <beta> and <gamma> are the priors.
            idx = l * K + k;
            var_beta[idx] = beta[idx] + zetabeta[idx] * var_beta_tmp[k];
            var_gamma[idx] = gamma[idx] + zetagamma[idx] * var_gamma_tmp[k];
        }
    }

    free( var_beta_tmp );
    free( var_gamma_tmp );
}

void P_update_logistic(const double* Dvarbeta, const double* Dvargamma, const double* mu, const double* Lambda, double* var_beta, double* var_gamma, double mintol, long L, long K)
{
    /*
    `Dvarbeta` and `Dvargamma` are a function of the values of `var_beta` and `var_gamma`. This
    dependence, however, is explicit on a set of latent populations assignments which in-turn
    depend on the variational parameters estimated at the previous step. So, the variables
    `Dvarbeta` and `Dvargamma` do not have to be updated.
    */

    long l, k, idx, numvar, update;
    long iter = 0;
    double tol = 10.0;
    double tmptol;
    double beta, gamma;
    double A, pbetagamma, pbeta, pgamma, ppbeta, ppgamma;
    double A_1, B_1, C_1, A_2, B_2, C_2;

    /*
    Iterate until succesive estimates are sufficiently similar
    or the number of iterations exceeds 1000.
    */
    while (tol>mintol && iter<1000) {

        numvar = 0;
        tol = 0.0;

        // loop over loci
        for (l=0; l<L; l++) {

            // only update a locus, if all its variational parameters 
            // satisfy positivity constraints.
            update = 1;
            for (k=0; k<K; k++) {
                if (var_beta[l*K+k]<=0 || var_gamma[l*K+k]<=0) {
                    update = 0;
                }
            }

            if (update==1) {

                // loop over populations
                for (k=0; k<K; k++) {

                    idx = l*K+k;
                    // compute pseudo-hyperparameters
                    pbetagamma = gsl_sf_psi_1(var_beta[idx]+var_gamma[idx]);
                    pbeta = gsl_sf_psi_1(var_beta[idx]);
                    pgamma = gsl_sf_psi_1(var_gamma[idx]);
                    ppbeta = gsl_sf_psi_n(2, var_beta[idx]);
                    ppgamma = gsl_sf_psi_n(2, var_gamma[idx]);

                    A_1 = pbeta-pbetagamma;
                    B_1 = -1.*pbetagamma;
                    A_2 = -1.*pbetagamma;
                    B_2 = pgamma-pbetagamma;
                    A = (gsl_sf_psi(var_beta[idx]) - gsl_sf_psi(var_gamma[idx]) - mu[l])*Lambda[k];
                    C_1 = -1.*A*pbeta - 0.5*Lambda[k]*ppbeta;
                    C_2 = A*pgamma - 0.5*Lambda[k]*ppgamma;

                    beta = (C_1*B_2-C_2*B_1)/(A_1*B_2-A_2*B_1);
                    gamma = (C_1*A_2-C_2*A_1)/(B_1*A_2-B_2*A_1);

                    // compute var_{beta,gamma}
                    tmptol = var_beta[idx];
                    var_beta[idx] = beta + Dvarbeta[idx];
                    tol += fabs(var_beta[idx]-tmptol);

                    tmptol = var_gamma[idx];
                    var_gamma[idx] = gamma + Dvargamma[idx];
                    tol += fabs(var_gamma[idx]-tmptol);

                    numvar += 1;

                }
            }
        }

        // compute convergence tolerance
        tol = 0.5*tol/numvar;
        iter += 1;
    }

    // printf("tol = %.8f, iter = %lu\n",tol,iter);
}
