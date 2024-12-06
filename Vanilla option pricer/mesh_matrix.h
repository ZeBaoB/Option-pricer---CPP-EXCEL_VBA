#pragma once
#include "pch.h"
#include "matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include <cmath>

/**
 * @file mesh_matrix.h
 * @brief Defines the `mesh_matrix` class, a specialized class for financial option pricing grids.
 *
 * @details
 * The `mesh_matrix` class extends the `matrix` class and represents a numerical grid used for option pricing.
 * It is parameterized by time and spot price mesh parameters, and it includes methods for solving pricing models
 * and retrieving option Greeks such as delta, gamma, theta, rho, and vega.
 *
 * @authors
 * - Bill Njonze
 * - Giorgio Enrico Maria Serva
 *
 * @date November 2024
 */

class mesh_matrix : public matrix {
private:
    matrix::matrix;
     long double T_;    ///< Maturity time of the option.
     long double s0_;   ///< Initial spot price of the underlying asset.
     long double K_;    ///< Strike price of the option.
     long double sinf_; ///< Maximum spot price (boundary condition).

public:
    /**
     * @brief Constructs a `mesh_matrix` with specified parameters.
     *
     * @param T Maturity time of the option.
     * @param s0 Initial spot price of the underlying asset.
     * @param K Strike price of the option.
     * @param spot_mest_parameter Number of spot price grid points (columns).
     * @param time_mesh_parameter Number of time grid points (rows).
     * @param sinf Maximum spot price (boundary condition).
     *
     * @details
     * Initializes a `mesh_matrix` as a `time_mesh_parameter x spot_mest_parameter` grid.
     * The parameters define the financial and numerical characteristics of the option pricing problem.
     */
    mesh_matrix( long double T,  long double s0,  long double K, unsigned int spot_mest_parameter, unsigned int time_mesh_parameter,  long double sinf)
        : matrix(time_mesh_parameter, spot_mest_parameter) {
        T_ = T;
        s0_ = s0;
        K_ = K;
        sinf_ = sinf;
    }

    /**
     * @brief Retrieves the maturity time of the option.
     *
     * @return The maturity time (`T_`).
     */
      long double  get_T() { return T_; }

    /**
     * @brief Retrieves the initial spot price of the underlying asset.
     *
     * @return The initial spot price (`s0_`).
     */
      long double  get_s0() { return s0_; }

    /**
     * @brief Retrieves the strike price of the option.
     *
     * @return The strike price (`K_`).
     */
      long double  get_K() { return K_; }

    /**
     * @brief Retrieves the maximum spot price (boundary condition).
     *
     * @return The maximum spot price (`sinf_`).
     */
    long double  get_sinf() { return sinf_; }

    /**
     * @brief Solves the option pricing model and populates the `mesh_matrix`.
     *
     * @param call Specifies whether the option is a call option (`true`) or a put option (`false`).
     * @param european Specifies whether the option is European (`true`) or American (`false`).
     * @param r_table A `matrix_plf` object representing the risk-free rate over time.
     * @param sigma The volatility of the underlying asset.
     */
    void solve(bool call, bool european, matrix_plf r_table,   long double  sigma);

    /**
     * @brief Retrieves the option value at a specific spot price.
     *
     * @param s0 The spot price at which to retrieve the option value.
     * @param line The time index (default is `1` for the initial time step).
     * @return The option value at the specified spot price and time index.
     */
    long double  retrieve_OptionValue( long double  s0, unsigned int line = 1);

    /**
     * @brief Retrieves the delta (sensitivity to spot price) of the option.
     *
     * @param s0 The spot price at which to compute the delta.
     * @return The delta value.
     */
    long double  retrieve_delta( long double s0);

    /**
     * @brief Retrieves the gamma (sensitivity of delta to spot price) of the option.
     *
     * @param s0 The spot price at which to compute the gamma.
     * @return The gamma value.
     */
    long double retrieve_gamma(  long double s0);

    /**
     * @brief Retrieves the theta (sensitivity to time) of the option.
     *
     * @param s0 The spot price at which to compute the theta.
     * @return The theta value (rate of change of option value with respect to time).
     */
      long double retrieve_theta(  long double s0);

    /**
     * @brief Retrieves the rho (sensitivity to interest rate) of the option.
     *
     * @param call Specifies whether the option is a call option (`true`) or a put option (`false`).
     * @param european Specifies whether the option is European (`true`) or American (`false`).
     * @param r_table A `matrix_plf` object representing the risk-free rate over time.
     * @param sigma The volatility of the underlying asset.
     * @param s0 The spot price at which to compute rho.
     * @return The rho value.
     */
      long double retrieve_rho(bool call, bool european, matrix_plf r_table,  long double sigma,   long double s0);

    /**
     * @brief Retrieves the vega (sensitivity to volatility) of the option.
     *
     * @param call Specifies whether the option is a call option (`true`) or a put option (`false`).
     * @param european Specifies whether the option is European (`true`) or American (`false`).
     * @param r_table A `matrix_plf` object representing the risk-free rate over time.
     * @param sigma The volatility of the underlying asset.
     * @param s0 The spot price at which to compute vega.
     * @return The vega value.
     */
      long double retrieve_vega(bool call, bool european, matrix_plf r_table,  long double sigma,   long double s0);
};
