#include "pch.h"
#include "matrix_tridiagonal.h"

matrix matrix_tridiagonal::inverse() const {
    unsigned int dim = get_number_of_lines();

    // Copy of the current matrix to work on for inversion
    matrix m = matrix(dim, dim);
    for (unsigned int i = 1; i <= dim; ++i) { for (unsigned int j = 1; j <= dim; ++j) { m(i, j) = (*this)(i, j); } }

    // Initialize the identity matrix to store the inverse
    matrix m_1(dim, dim);
    for (unsigned int i = 1; i <= dim; ++i) {
        m_1(i, i) = 1.0;
    }

    // Lower triangular elimination
    for (unsigned int i = 1; i < dim; ++i) {
          long double  mii = m(i, i);
          long double  alpha = m(i + 1, i);

        m.add_line_to_line(i + 1, -alpha / mii, i);
        m_1.add_line_to_line(i + 1, -alpha / mii, i);
    }

    // Upper triangular elimination
    for (unsigned int i = dim; i > 1; --i) {
          long double  mii = m(i, i);
          long double  alpha = m(i - 1, i);

        m.add_line_to_line(i - 1, -alpha / mii, i);
        m_1.add_line_to_line(i - 1, -alpha / mii, i);
    }

    // Normalize diagonal to 1
    for (unsigned int i = 1; i <= dim; ++i) {
          long double  mii = m(i, i);
        for (unsigned int j = 1; j <= dim; ++j) {
            m_1(i, j) /= mii;
        }
    }
    return m_1;
};

matrix matrix_tridiagonal::gauss_saidel(matrix b, matrix x0) const {

    unsigned int dim = b.get_number_of_lines();
    matrix x = x0;      // Initialize x0

    // Repetition for convergence
    const unsigned int max_iter = 1000; // Max number of iteration 
    const   long double  tol = 1e-12;   // Tollerance for convergence 
    for (unsigned int iter = 0; iter < max_iter; ++iter) {
          long double  max_diff = 0.0; // Difference between old and new x

        for (unsigned int i = 1; i <= dim; ++i) {
              long double  sum = b(i, 1); // Initialize the term

            // Diagonal element
              long double  a_ii = (*this)(i, i);

            // Sum of elements on the left
            if (i > 1) {
                  long double  a_ij = (*this)(i, i - 1);
                sum -= a_ij * x(i - 1, 1);
            }

            // Sum of the elements on the right
            if (i < dim) {
                  long double  a_ij = (*this)(i, i + 1);
                sum -= a_ij * x0(i + 1, 1);
            }

            // Update the solution
              long double  x_new = sum / a_ii;
            max_diff = (max_diff > std::abs(x_new - x(i, 1))) ? max_diff : std::abs(x_new - x(i, 1));
            x(i, 1) = x_new;
        }

        // Control the convergence
        if (max_diff < tol) {
            break;
        }

        // Update x0 for the next iteration
        x0 = x;
    }
    return x;
}