#include "Matrix.h"
#include <iostream>
#include <stdexcept>

int main() {
    /*
    matrix m1(5, 6);
    for (int i = 1; i <= 5; i++) {
        for (int j = i; j <= 6; j++) {
            m1(i, j) = 1.0;
        }
    }

    matrix m2(5, 6);
    for (int i = 1; i <= 5; i++) {
        for (int j = 1; j <= i; j++) {
            m2(i, j) = 1.0;
        }
    }

    matrix m3(6, 5);
    m3(1, 1) = 1.0;
    m3(2, 2) = 2.0;
    m3(3, 3) = 3.0;
    m3(4, 4) = 4.0;
    m3(5, 5) = 5.0;

    std::cout << m1 << std::endl << std::endl;
    std::cout << m2 << std::endl << std::endl;
    std::cout << m1 * m3 << std::endl;
    std::cout << m2 * m3 << std::endl;


    std::cout << returnTable << std::endl << std::endl;

    //Example lookup
    try {
        double t = 1.5;
        std::cout << "Value at t = " << t << ": " << returnTable(t) << std::endl;
    }
    catch (const std::out_of_range& e) {
        std::cerr << e.what() << std::endl;
    }

    double diag[4]{ 3, 3, 3, 3 };
    double diagSup[3]{ -1, -1, -1 };
    double diagInf[3]{ -1, -1, -1 };

    matrix_tridiagonal m_tri(4, diag, diagSup, diagInf);
    matrix m_tri_inv = m_tri.inverse();

    std::cout << m_tri << std::endl << std::endl;
    std::cout << m_tri_inv << std::endl << std::endl;
    std::cout << m_tri * m_tri.inverse() << std::endl << std::endl;
    std::cout << m_tri.inverse() * m_tri << std::endl << std::endl;
    */
    bool call = false;
    bool european = true;
    double maturity = 1.5;
    double s0 = 42.3;
    double k = 45.0;
    double tab_t[4]{ 0.0, 1.0, 1.75, 2 };
    double tab_r[4]{ 0.025, 0.025, 0.025, 0.030 };
    double sigma = 0.25;
    int spot_mesh_parameter = 50;
    int time_mesh_parameter = 50;

    double sinf = s0 *5; // 15 must depend of spot_mesh_parameter 

    matrix_plf returnTable = matrix_plf(tab_t, tab_r, 4);
    //std::cout << returnTable(50.2) << std::endl;

    //std::cout << "Boundary condition at S0=0: " << k * std::exp(-returnTable.integral(0.0, maturity)) << std::endl << std::endl;

    
    mesh_matrix m_matrix(maturity, s0, k, spot_mesh_parameter, time_mesh_parameter, sinf);
    
    //std::cout << m_matrix << std::endl << std::endl;
    m_matrix.solve(call, european, returnTable, sigma);
    for (int j = 1; j <= 5; ++j) {
        std::cout << "Value at S = "<< sinf * (j-1)/ spot_mesh_parameter << ": " << m_matrix(1, j) << std::endl;
    }
    //std::cout << m_matrix << std::endl << std::endl;
    std::cout << "Price : " << m_matrix.retrieve_OptionValue(s0) << std::endl;
    std::cout << "Delta : " << m_matrix.retrieve_delta(s0) << std::endl;
    std::cout << "Gamma : " << m_matrix.retrieve_gamma(s0) << std::endl;
    std::cout << "Theta : " << m_matrix.retrieve_theta(s0) << std::endl;

    std::cout << "Rho : " << m_matrix.retrieve_rho(call, european, returnTable, sigma, s0) << std::endl;
    std::cout << "vega : " << m_matrix.retrieve_vega(call, european, returnTable, sigma, s0) << std::endl;
    

    return 0;
}