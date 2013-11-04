
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

        // compute digamma functions of parameters
        for (k=0; k<K; k++) {

            // initialize var_{beta,gamma}
            var_beta_tmp[k] = 0.0;
            var_gamma_tmp[k] = 0.0;
        }

        // sum over samples
        for (n=0; n<N; n++) {

            genotype = G[n*L+l];
            if (genotype!=3) {

                // compute xi*theta_{beta,gamma}
                theta_beta_sum = 0.0;
                theta_gamma_sum = 0.0;
                for (k=0; k<K; k++) {
                    theta_beta_sum += xi[n*K+k] * zetabeta[l*K+k];
                    theta_gamma_sum += xi[n*K+k] * zetagamma[l*K+k];
                }

                // increment var_{beta,gamma}
                for (k=0; k<K; k++) {
                    var_beta_tmp[k] += (double) genotype * xi[n*K+k] / theta_beta_sum;
                    var_gamma_tmp[k] += (double) (2-genotype) * xi[n*K+k] / theta_gamma_sum;
                }
            }
        }

        // compute var_{beta,gamma}
        for (k=0; k<K; k++) {
            // change `hyper` for F model
            idx = l*K+k;
            var_beta[idx] = beta[idx] + zetabeta[idx] * var_beta_tmp[k];
            var_gamma[idx] = gamma[idx] + zetagamma[idx] * var_gamma_tmp[k];
        }
    }

    free( var_beta_tmp );
    free( var_gamma_tmp );
}

void P_update_logistic(const double* Dvarbeta, const double* Dvargamma, const double* mu, const double* Lambda, double* var_beta, double* var_gamma, double mintol, long L, long K)
{
    // ZETABETA and ZETAGAMMA should be fixed since it is actually a component of latent variable Z that is kept fixed in this step.

    long l, k, idx, numvar, update;
    long iter = 0;
    double tol = 10.0;
    double tmptol;
    double beta, gamma;
    double A, pbetagamma, pbeta, pgamma, ppbeta, ppgamma;
    double A_1, B_1, C_1, A_2, B_2, C_2;

    while (tol>mintol && iter<1000) {

        numvar = 0;
        tol = 0.0;
        // loop over loci
        for (l=0; l<L; l++) {

            update = 1;
            for (k=0; k<K; k++) {
                if (var_beta[l*K+k]<=0 || var_gamma[l*K+k]<=0) {
                    update = 0;
                }
            }

            // only update a locus, if all its variational parameters are not at constraint boundaries.
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

        tol = 0.5*tol/numvar;
        iter += 1;
    }

    // printf("tol = %.8f, iter = %lu\n",tol,iter);
}
