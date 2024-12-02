#include "pch.h"
#include "pricer.h"
#include "matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include "mesh_matrix.h"
double computePriceAndPartialDerivates(
    double* result,
    double* underlying,
    double* option_price,
    double* delta,
    double* exerciceBoundary,
    bool isCall,
    bool isEuropean,
    double maturity,
    double strike,
    int time_mesh_params,
    int spot_mesh_params,
    double S0,
    const double* time_steps,  // Time steps array
    const double* r_values,   // Corresponding interest rates array
    int size,                 // Size of the r table (number of time steps)
    double sigma)
{

    double sinf = S0 * 2.5;
    long double ds = sinf / (spot_mesh_params - 1);

    matrix_plf returnTable = matrix_plf(time_steps, r_values, size);
    mesh_matrix m_matrix(maturity, S0, strike, spot_mesh_params, time_mesh_params, sinf);

    m_matrix.solve(isCall, isEuropean, returnTable, sigma);

    // Initialize arrays
    underlying[0] = 0.0;
    option_price[0] = m_matrix(1, 1);

    for (int j = 2; j < spot_mesh_params - 2; ++j) {
        underlying[j - 1] = (j - 1) * ds;
        delta[j - 1] = m_matrix.retrieve_delta((j - 1) * ds);
        option_price[j - 1] = m_matrix(1, j);
    }

    for (int i = 1; i <= time_mesh_params - 2; ++i) {
        exerciceBoundary[i - 1] = isCall ? sinf : 0.0;

        for (int jj = isCall ? 1 : spot_mesh_params;
            isCall ? (jj <= spot_mesh_params) : (jj >= 1);
            jj += isCall ? 1 : -1) {
            double S = sinf * (jj - 1) / (spot_mesh_params - 1);
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
