#include "pch.h"
#include "Matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include "mesh_matrix.h"

#include <iostream>
#include <stdexcept>

#include <chrono> // Pour mesurer le temps d'exécution

/**
 * @file Test.cpp
 * @brief Code used to test the methods created and the accuracy of the algorithm
 *
 * @details
 * In order to validate the results we compare the numerical method with the explicit solution given by Black–Scholes formula in the case where r is constant.
 * We also verify the equivalence between american and european exercises for calls with r > 0 and for puts with r < 0.
 *
 * @authors
 * - Bill Njonze
 * - Giorgio Enrico Maria Serva
 *
 * @date November 2024
 */

int main() {

    bool call = true;
    bool european = true;
    double maturity = 0.5;
    double s0 = 42.0;
    double k = 45.0;
    double tab_t[4]{ -0.1, 0.2, 0.3, 0.4 };
    double tab_r[4]{ 0.05, 0.05, 0.05, 0.05 };
    double sigma = 0.25;
    int spot_mesh_parameter = 250;
    int time_mesh_parameter = 250;

    double sinf = s0 * 3;

    matrix_plf returnTable = matrix_plf(tab_t, tab_r, 4);

    std::cout << "----- Constant risk free rate ----- Spot = 42, Strike = 45, Maturity = 6 mois, r = 0.5, sigma = 25% -----------------------------------------------------\n\n";

    std::cout << "European call option (spot_mesh_parameter = 250, time_mesh_parameter = 250) : \n";
    // Début de la mesure du temps
    auto start = std::chrono::high_resolution_clock::now();
    mesh_matrix m_matrix(maturity, s0, k, spot_mesh_parameter, time_mesh_parameter, sinf);
    m_matrix.solve(call, european, returnTable, sigma);
    std::cout << "Computed Price : " << m_matrix.retrieve_OptionValue(s0) << " VS 2.17311 (the explicit solution)" << std::endl;
    std::cout << "Computed Delta : " << m_matrix.retrieve_delta(s0) << " VS 0.43625 (the explicit solution)" << std::endl;
    std::cout << "Computed Gamma : " << m_matrix.retrieve_gamma(s0) << " VS 0.05304 (the explicit solution)" << std::endl;
    std::cout << "Computed Theta : " << m_matrix.retrieve_theta(s0) << " VS 3.73158 (the explicit solution)" << std::endl;
    std::cout << "Computed Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << " VS 8.074791 (the explicit solution)" << std::endl;
    std::cout << "Computed Vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << " VS 11.69641 (the explicit solution)" << std::endl;
    // Fin de la mesure du temps
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl << std::endl << std::endl;


    std::cout << "European put option (spot_mesh_parameter = 250, time_mesh_parameter = 250) : \n";
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    m_matrix = mesh_matrix(maturity, s0, k, spot_mesh_parameter, time_mesh_parameter, sinf);
    m_matrix.solve(false, european, returnTable, sigma);
    std::cout << "Computed Price : " << m_matrix.retrieve_OptionValue(s0) << " VS 4.06205 (the explicit solution)" << std::endl;
    std::cout << "Computed Delta : " << m_matrix.retrieve_delta(s0) << " VS -0.56375 (the explicit solution)" << std::endl;
    std::cout << "Computed Gamma : " << m_matrix.retrieve_gamma(s0) << " VS 0.05304 (the explicit solution)" << std::endl;
    std::cout << "Computed Theta : " << m_matrix.retrieve_theta(s0) << " VS 1.53713 (the explicit solution)" << std::endl;
    std::cout << "Computed Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << " VS -13.86969 (the explicit solution)" << std::endl;
    std::cout << "Computed Vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << " VS 11.69641 (the explicit solution)" << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
    std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl << std::endl;


    returnTable.matrix::operator()(2, 3) = 0.04;
    returnTable.matrix::operator()(2, 4) = 0.04;
    std::cout << "risk free rate at time = 0.25 : " << returnTable(0.25) << std::endl;
    std::cout << "risk free rate at time = 0.30 : " << returnTable(0.30) << std::endl;
    std::cout << "-------------------------- Call ------------------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "-------------------------- Spot = 42, Strike = 45, Maturity = 6 mois, r =environ 0.45 (positive), sigma = 25% ------------------------------------------------------\n\n";

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "European. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(call, european, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;


    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "\nEuropean. spot_mesh_parameter = 300, time_mesh_parameter = 300 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 300, 300, sinf);
    m_matrix.solve(call, european, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;


    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "\nAmerican. spot_mesh_parameter = 300, time_mesh_parameter = 300 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 300, 300, sinf);
    m_matrix.solve(call, false, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;



    std::cout << "\n\n\n\n-------------------------- Put -------------------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "-------------------------- Spot = 42, Strike = 45, Maturity = 6 mois, r =environ 0.45 (positive), sigma = 25% ------------------------------------------------------\n\n";

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "European. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(false, true, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "\nAmerican. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(false, false, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;




    //Negative r
    std::cout << "\n\n\n\n-------------------------- Negative risk free interest rate  -----------------------------------------------------------------------------------------------------\n\n";
    for (int i = 1; i <= 4; ++i) { returnTable.matrix::operator()(2, i) *= -1; }
    std::cout << "\n\nrisk free rate at time = 0.025 : " << returnTable(0.25) << std::endl << std::endl;


    std::cout << "-------------------------- Call -------------------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "-------------------------- Spot = 42, Strike = 45, Maturity = 6 mois, r =environ -0.45 (negative), sigma = 25% ------------------------------------------------------\n\n";

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "\nEuropean. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(call, european, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "\nAmerican. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(call, false, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;



    std::cout << "\n\n\n\n-------------------------- Put -------------------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "-------------------------- Spot = 42, Strike = 45, Maturity = 6 mois, r =environ -0.45 (negative), sigma = 25% ------------------------------------------------------\n\n";

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "European. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(false, true, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    //---------------------------------------------------------------------------------------------------------------
    // Début de la mesure du temps
    start = std::chrono::high_resolution_clock::now();
    std::cout << "\nAmerican. spot_mesh_parameter = 100, time_mesh_parameter = 100 \n";
    m_matrix = mesh_matrix(maturity, s0, k, 100, 100, sinf);
    m_matrix.solve(false, false, returnTable, sigma);
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;
    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    // Fin de la mesure du temps
    end = std::chrono::high_resolution_clock::now(); elapsed = end - start;
    // Affichage du temps d'exécution
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}