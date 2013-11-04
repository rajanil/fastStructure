
typedef unsigned char uint8_t;

void P_update_simple( const uint8_t* G, const double* zetabeta, const double* zetagamma, const double* xi, const double* beta, const double* gamma, double* var_beta, double* var_gamma, long N, long L, long K );

void P_update_logistic( const double* Dvarbeta, const double* Dvargamma, const double* mu, const double* Lambda, double* var_beta, double* var_gamma, double mintol, long L, long K);
