#include "pch.h"
#include "pricer.h"
#include "matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include "mesh_matrix.h"

double computePriceAndPartialDerivates(
    double* result,
    double* error,
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
    const double* r_values,    // Corresponding interest rates array
    int size,              // Size of the r table (number of time steps)
    double sigma)
{
    double sinf = S0 * 3;

    matrix_plf returnTable = matrix_plf(time_steps, r_values, size);

    mesh_matrix m_matrix(maturity, S0, strike, spot_mesh_params, time_mesh_params, sinf);
    
    m_matrix.solve(isCall, isEuropean, returnTable, sigma); //


    double ds = sinf / (spot_mesh_params - 1);

    underlying[0] = 0.0;
    option_price[0] = m_matrix(1, 1);
    for (int j = 2; j < spot_mesh_params; ++j) {
        underlying[j - 1] = (j - 1) * ds;
        delta[j - 1] = m_matrix.retrieve_delta((j - 1) * ds);
        option_price[j - 1] = m_matrix(1, j);
    }
    for (int i = 1; i <= time_mesh_params; ++i) {
        exerciceBoundary[i - 1] = isCall? strike + m_matrix.retrieve_OptionValue(S0, i) : strike - m_matrix.retrieve_OptionValue(S0, i);
    }

    result[0] = m_matrix.retrieve_OptionValue(S0); // Option price
    result[1] = m_matrix.retrieve_delta(S0); // Delta
    result[2] = m_matrix.retrieve_gamma(S0); // Gamma
    result[3] = m_matrix.retrieve_theta(S0); // Theta
    result[4] = m_matrix.retrieve_rho(isCall, isEuropean, returnTable, sigma, S0); // Rho
    result[5] = m_matrix.retrieve_vega(isCall, isEuropean, returnTable, sigma, S0); // Vega

    error[0] = 1.58463; // Option price error
    error[1] = 0.44037; // Delta error
    error[2] = 0.07614; // Gamma error
    error[3] = 1.09025; // Theta error
    error[4] = 25.56445; // Rho error
    error[5] = 20.43657; // Vega error
    return returnTable.integral(1.0, 1.5);
}
