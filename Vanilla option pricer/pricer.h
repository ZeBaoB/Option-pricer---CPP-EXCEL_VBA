#ifdef VanillaOptionPricer_EXPORTS
#define VanillaOptionPricer __declspec(dllexport)
#else
#define VanillaOptionPricer __declspec(dllimport)
#endif

extern "C" VanillaOptionPricer double computePriceAndPartialDerivates(
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
;