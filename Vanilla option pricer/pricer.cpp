#include "pch.h"
#include "pricer.h"
#include "matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include "mesh_matrix.h"
long double computePriceAndPartialDerivates(
    long double* result,
    long double* underlying,
    long double* option_price,
    long double* delta,
    long double* exerciceBoundary,
    bool isCall,
    bool isEuropean,
    long double maturity,
    long double strike,
    unsigned int time_mesh_params,
    unsigned int spot_mesh_params,
    long double S0,
    const long double* time_steps,  // Time steps array
    const long double* r_values,   // Corresponding interest rates array
    unsigned int size,             // Size of the r table (number of time steps)
    long double sigma)
{
    long double etalon = (S0 > strike) ? S0 : strike;
    long double sinf = 1.5*etalon + 3 * S0 * sigma * maturity;
    long double ds = sinf / (spot_mesh_params - 1);

    matrix_plf returnTable = matrix_plf(time_steps, r_values, size);
    mesh_matrix m_matrix(maturity, S0, strike, spot_mesh_params, time_mesh_params, sinf);

    m_matrix.solve(isCall, isEuropean, returnTable, sigma);

    // Initialize arrays
    underlying[0] = 0.0;
    option_price[0] = m_matrix(1, 1);

    for (unsigned int j = 2; j < spot_mesh_params - 2; ++j) {
        underlying[j - 1] = (j - 1) * ds;
        delta[j - 1] = m_matrix.retrieve_delta((j - 1) * ds);
        option_price[j - 1] = m_matrix(1, j);
    }

    for (unsigned int i = 1; i <= time_mesh_params - 2; ++i) {
        exerciceBoundary[i - 1] = isCall ? sinf : 0.0;

        for (unsigned int jj = isCall ? 1 : spot_mesh_params;
            isCall ? (jj <= spot_mesh_params) : (jj >= 1);
            jj += isCall ? 1 : -1) {
            long double  S = sinf * (jj - 1) / (spot_mesh_params - 1);
            if ((isCall ? (S - strike) : (strike - S)) - m_matrix(i, jj) == 0) {
                exerciceBoundary[i - 1] = S;
                break;
            }
        }
    }
    result[0] = m_matrix.retrieve_OptionValue(S0);
    result[1] = m_matrix.retrieve_delta(S0);
    result[2] = m_matrix.retrieve_gamma(S0);
    result[3] = m_matrix.retrieve_theta(S0);
    result[4] = m_matrix.retrieve_rho(isCall, isEuropean, returnTable, sigma, S0);
    result[5] = m_matrix.retrieve_vega(isCall, isEuropean, returnTable, sigma, S0);

    return 0.0;
}