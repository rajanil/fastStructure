
#include "marglikehood.h"
#include <math.h>

double marglikehood(const uint8_t* G, const double* zetabeta, const double* zetagamma, const double* xi, long N, long L, long K)
{
    uint8_t genotype;
    long n, l, k;
    double E, zasum, zbsum;

    E = 0;
    // loop over loci
    for (l=0; l<L; l++) {

        // loop over samples 
        for (n=0; n<N; n++) {

            genotype = G[n*L+l];

            // missing data do not contribute
            if (genotype!=3) {

                zasum = 0.;
                zbsum = 0.;

                // loop over populations
                if (genotype==0) {

                    for (k=0; k<K; k++) {
                        zasum += zetagamma[l*K+k]*xi[n*K+k];
                    }
                    zbsum = zasum;

                } else if (genotype==1) {

                    for (k=0; k<K; k++) {
                        zasum += zetagamma[l*K+k]*xi[n*K+k];
                        zbsum += zetabeta[l*K+k]*xi[n*K+k];
                    }

                } else if (genotype==2) {

                    for (k=0; k<K; k++) {
                        zasum += zetabeta[l*K+k]*xi[n*K+k];
                    }
                    zbsum = zasum;

                }

                E += log(zasum) + log(zbsum);
            }

        }
    }
    return E;
}
