#ifdef VanillaOptionPricer_EXPORTS
#define VanillaOptionPricer __declspec(dllexport)
#else
#define VanillaOptionPricer __declspec(dllimport)
#endif

/**
 * @file pricer.h
 * @brief Exports the `computePriceAndPartialDerivates` function for option pricing and Greeks computation.
 *
 * @details
 * This header defines the `computePriceAndPartialDerivates` function as an exported or imported symbol
 * depending on the build configuration. The function is compatible with C linkage for easier integration
 * with other programming languages and systems.
 */

 /**
  * @fn VanillaOptionPricer double computePriceAndPartialDerivates(
  *     double* result,
  *     double* underlying,
  *     double* option_price,
  *     double* delta,
  *     double* exerciceBoundary,
  *     bool isCall,
  *     bool isEuropean,
  *     double maturity,
  *     double strike,
  *     int time_mesh_params,
  *     int spot_mesh_params,
  *     double S0,
  *     const double* time_steps,
  *     const double* r_values,
  *     int size,
  *     double sigma)
  * 
 * @brief Defines the `computePriceAndPartialDerivates` function for pricing options and computing their partial derivatives.
 *
 * @details
 * This function computes the price of an option and its Greeks (delta, gamma, theta, rho, and vega)
 * based on the specified parameters. It constructs a mesh grid using `mesh_matrix` to solve the
 * option pricing model and retrieves the necessary values.
 *
 * @param[out] result Array to store computed option price and partial derivatives (size 6):
 *  - `result[0]`: Option price.
 *  - `result[1]`: Delta.
 *  - `result[2]`: Gamma.
 *  - `result[3]`: Theta.
 *  - `result[4]`: Rho.
 *  - `result[5]`: Vega.
 * @param[out] underlying Array to store the underlying spot price values along the grid.
 * @param[out] option_price Array to store the computed option prices along the grid.
 * @param[out] delta Array to store the delta values along the grid.
 * @param[out] exerciceBoundary Array to store the exercise boundary at different time steps.
 * @param[in] isCall Specifies if the option is a call option (`true`) or put option (`false`).
 * @param[in] isEuropean Specifies if the option is European (`true`) or American (`false`).
 * @param[in] maturity The maturity time of the option.
 * @param[in] strike The strike price of the option.
 * @param[in] time_mesh_params Number of time grid points.
 * @param[in] spot_mesh_params Number of spot price grid points.
 * @param[in] S0 The initial spot price of the underlying asset.
 * @param[in] time_steps Array of time steps for the risk-free rate table.
 * @param[in] r_values Array of corresponding risk-free interest rates.
 * @param[in] size The size of the `time_steps` and `r_values` arrays.
 * @param[in] sigma The volatility of the underlying asset.
 *
 * @return Returns `0.0` upon successful computation.
 *
 * @details
 * 1. Constructs a `matrix_plf` object to represent the risk-free rate over time.
 * 2. Constructs a `mesh_matrix` object to model the option pricing grid.
 * 3. Uses the `mesh_matrix` object to solve the option pricing problem.
 * 4. Retrieves and stores the option price and Greeks in the `result` array.
 * 5. Populates the `underlying`, `option_price`, `delta`, and `exerciceBoundary` arrays.
 *
 * @note
 * - The spot price and time grids are evenly spaced based on the parameters provided.
 * - The Greeks are computed at the specified spot price `S0` and along the spot price grid.
 */
extern "C" VanillaOptionPricer double computePriceAndPartialDerivates(
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
    const double* r_values,    // Corresponding interest rates array
    int size,              // Size of the r table (number of time steps)
    double sigma)
;