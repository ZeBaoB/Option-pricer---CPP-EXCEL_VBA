#include "pch.h"
#include "matrix_tridiagonal.h"

matrix matrix_tridiagonal::inverse() const {
    int dim = get_number_of_lines();

    // Copy of the current matrix to work on for inversion
    matrix m = matrix(dim, dim);
    for (int i = 1; i <= dim; ++i) { for (int j = 1; j <= dim; ++j) { m(i, j) = (*this)(i, j); } }

    // Initialize the identity matrix to store the inverse
    matrix m_1(dim, dim);
    for (int i = 1; i <= dim; ++i) {
        m_1(i, i) = 1.0;
    }

    // Lower triangular elimination
    for (int i = 1; i < dim; ++i) {
        double mii = m(i, i);
        double alpha = m(i + 1, i);

        m.add_line_to_line(i + 1, -alpha / mii, i);
        m_1.add_line_to_line(i + 1, -alpha / mii, i);
    }

    // Upper triangular elimination
    for (int i = dim; i > 1; --i) {
        double mii = m(i, i);
        double alpha = m(i - 1, i);

        m.add_line_to_line(i - 1, -alpha / mii, i);
        m_1.add_line_to_line(i - 1, -alpha / mii, i);
    }

    // Normalize diagonal to 1
    for (int i = 1; i <= dim; ++i) {
        double mii = m(i, i);
        for (int j = 1; j <= dim; ++j) {
            m_1(i, j) /= mii;
        }
    }
    return m_1;
};