#include "Matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include "mesh_matrix.h"

#include <iostream>
#include <stdexcept>

#include <chrono> // Pour mesurer le temps d'exécution

int main() {
    // Début de la mesure du temps
    auto start = std::chrono::high_resolution_clock::now();

    bool call = false;
    bool european = false;
    double maturity = 1.5;
    double s0 = 42.3;
    double k = 45.0;
    double tab_t[4]{ 0.0, 1.0, 1.75, 2 };
    double tab_r[4]{ -0.025, -0.025, -0.025, -0.025 };
    double sigma = 0.25;
    int spot_mesh_parameter = 500;
    int time_mesh_parameter = 500;

    double sinf = s0 * 3; 

    matrix_plf returnTable = matrix_plf(tab_t, tab_r, 4);

    mesh_matrix m_matrix(maturity, s0, k, spot_mesh_parameter, time_mesh_parameter, sinf);
    // std::cout << m_matrix << std::endl << std::endl;
    m_matrix.solve(call, european, returnTable, sigma);


    std::cout << "spot_mesh_parameter : " << spot_mesh_parameter << std::endl;
    std::cout << "time_mesh_parameter : " << time_mesh_parameter << std::endl << std::endl;

    //std::cout << m_matrix << std::endl << std::endl;
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;

    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    

    // Fin de la mesure du temps
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}