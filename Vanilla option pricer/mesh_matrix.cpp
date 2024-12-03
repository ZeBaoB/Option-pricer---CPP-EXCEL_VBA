#include "pch.h"
#include "mesh_matrix.h"
#include <cmath>

void mesh_matrix::solve(bool call, bool european, matrix_plf r_table, double sigma) {
    int dim = get_number_of_columns();
    int nl = get_number_of_lines();
    matrix Vi(dim, 1);
    matrix Vi_1(dim, 1);
    long double r;
    long double ak;
    long double bk;
    long double ck;
    long double sk;
    long double dt = T_ / (nl - 1);
    long double S_K; // S-K or K-S
    double coef = call ? 1.0 : -1.0;

    //Set the boundary conditions of the mesh
    if (call) {
        // If the spot price is 0, the price of the call is zero
        fill_column(1, 0.0);
        // If the spot price is Sinfini, the price of the call is S - K * exp(-rt)
        for (int i = 1; i <= nl; i++) {
            matrix::operator()(i, dim) = get_sinf() - get_K() * std::exp(-r_table.integral(get_T() * (i - 1) / (nl - 1), get_T()));
        }
        // At time t = T, the price is the payoff
        for (int j = 2; j <= dim; j++) {
            S_K = get_sinf() * (j - 1) / (dim - 1) - get_K();
            matrix::operator()(nl, j) = (S_K < 0.0) ? 0.0 : S_K;
        }
    }
    else {
        // If the spot price is Sinfini, the price of the put is zero
        fill_column(dim, 0.0);
        // If the spot price is 0, the price of the put is K * exp(-rt)
        for (int i = 1; i <= nl; i++) {
            matrix::operator()(i, 1) = get_K() * std::exp(-r_table.integral(get_T() * (i - 1) / (nl - 1), get_T()));
        }
        // At time t = T, the price is the payoff
        for (int j = 2; j <= dim; j++) {
            S_K = get_K() - get_sinf() * (j - 1) / (dim - 1);
            matrix::operator()(nl, j) = (S_K < 0.0) ? 0.0 : S_K;
        }
    }

    // The rows represent time from t0 to T. We go backward to retrieve the price
    for (int i = nl; i >= 2; --i) {
        for (int k = 1; k <= dim; ++k) {
            Vi(k, 1) = matrix::operator()(i, k);
        }

        r = (r_table(T_ * (i - 1) / nl) + r_table(T_ * i / nl)) / 2;

        matrix_tridiagonal Lmatrix(dim);
        matrix_tridiagonal Rmatrix(dim);

        Lmatrix(1, 1) = (matrix::operator()(i - 1, 1) > 0) ? 1 / matrix::operator()(i - 1, 1) : 1;
        Rmatrix(1, 1) = (Vi(1, 1) > 0.0) ? 1.0 / Vi(1, 1) : 1;
        Lmatrix(dim, dim) = (matrix::operator()(i - 1, dim) > 0) ? 1 / matrix::operator()(i - 1, dim) : 1;
        Rmatrix(dim, dim) = (Vi(dim, 1) > 0.0) ? 1.0 / Vi(dim, 1) : 1;

        //Compute the coefficients of matrix
        for (int k = 2; k <= dim - 1; k++) {
            sk = 1.0 * (k - 1);
            ak = dt / 4 * (sigma * sigma * sk * sk - r * sk);
            bk = -dt / 2 * (sigma * sigma * sk * sk + r);
            ck = dt / 4 * (sigma * sigma * sk * sk + r * sk);
            Lmatrix(k, k - 1) = -ak;
            Lmatrix(k, k) = 1 - bk;
            Lmatrix(k, k + 1) = -ck;
            Rmatrix(k, k - 1) = ak;
            Rmatrix(k, k) = 1 + bk;
            Rmatrix(k, k + 1) = ck;
        }
        //Vi_1 = Lmatrix.inverse() * Rmatrix * Vi;
        Vi_1 = Lmatrix.gauss_saidel(Rmatrix * Vi, Vi);

        if (european) { // european call or put option
            for (int k = 1; k <= dim; ++k) {
                matrix::operator()(i - 1, k) = Vi_1(k, 1);
            }
        }
        else { // american
            for (int k = 1; k <= dim ; ++k) {
                matrix::operator()(i - 1, k) = (Vi_1(k, 1) > coef * (get_sinf() * (k - 1) / (dim - 1) - get_K())) ? Vi_1(k, 1) : coef * (get_sinf() * (k - 1) / (dim - 1) - get_K());
            }
        }
    }
}

long double mesh_matrix::retrieve_OptionValue(long double s0, int line)
{
    double j;
    int j1;
    int j2;
    j = 1.0 + s0 * (get_number_of_columns() - 1) / get_sinf();
    j1 = static_cast<int>(std::floor(j));
    j2 = static_cast<int>(std::floor(j) + 1);

    return matrix::operator()(line, j1) + (j - j1) * (matrix::operator()(line, j2) - matrix::operator()(line, j1));
}

long double mesh_matrix::retrieve_delta(long double s0)
{
    long double ds = get_sinf() / (get_number_of_columns() - 1);
    return (retrieve_OptionValue(s0 + ds) - retrieve_OptionValue(s0 - ds)) / 2 / ds;
}

long double mesh_matrix::retrieve_gamma(long double s0)
{
    long double ds = get_sinf() / (get_number_of_columns() - 1);
    return (retrieve_OptionValue(s0 + ds) + retrieve_OptionValue(s0 - ds) - 2 * retrieve_OptionValue(s0)) / ds / ds;
}

long double mesh_matrix::retrieve_theta(long double s0)
{
    return ( retrieve_OptionValue(s0, 2) - retrieve_OptionValue(s0) ) / get_T() * (get_number_of_lines() - 1);
}

long double mesh_matrix::retrieve_rho(bool call, bool european, matrix_plf r_table, double sigma, long double s0)
{
    long double price1 = retrieve_OptionValue(s0);
    matrix_plf new_rtable = r_table.shift_value(1e-4);
    solve(call, european, new_rtable, sigma);
    long double price2 = retrieve_OptionValue(s0);

    solve(call, european, r_table, sigma);

    return (price2 - price1) / 1e-4;
}

long double mesh_matrix::retrieve_vega(bool call, bool european, matrix_plf r_table, double sigma, long double s0)
{
    long double price1 = retrieve_OptionValue(s0);
    solve(call, european, r_table, sigma + 1e-4);
    long double price2 = retrieve_OptionValue(s0);

    solve(call, european, r_table, sigma);

    return (price2 - price1) / 1e-4;
}

