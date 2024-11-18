#include "pch.h"
#include "matrix_plf.h"
#include <stdexcept>


double matrix_plf::operator()(double t) const {
    // Ensure `t` is within the bounds of `tab_t` for interpolation
    int columns = get_number_of_columns();
    for (int j = 1; j < columns; ++j) {
        double t1 = matrix::operator()(1, j);
        double t2 = matrix::operator()(1, j + 1);
        if (t >= t1 && t <= t2) {
            double r1 = matrix::operator()(2, j);
            double r2 = matrix::operator()(2, j + 1);
            // Linear interpolation
            return r1 + (r2 - r1) * (t - t1) / (t2 - t1);
        }
    }
    // Throw exception if `t` is out of bounds
    throw std::out_of_range("Value of t is out of bounds.");
}


matrix_plf matrix_plf::shift_value(double delta_r)
{
    int numberOfPoints = get_number_of_columns();

    // Create arrays to hold current matrix values
    double* tab_t = new double[numberOfPoints];
    double* tab_r = new double[numberOfPoints];

    // Fill `tab_t` and `tab_r` with the current values
    for (int j = 1; j <= numberOfPoints; j++) {
        tab_t[j - 1] = matrix::operator()(1, j);          // Copy the first row (t values)
        tab_r[j - 1] = matrix::operator()(2, j) + delta_r; // Shift the second row (r values) by `delta_r`
    }

    // Create a new `matrix_plf` with the shifted values
    matrix_plf m(tab_t, tab_r, numberOfPoints);

    // Clean up dynamically allocated memory
    delete[] tab_t;
    delete[] tab_r;

    return m;
}


double matrix_plf::integral(double t, double T) const
{
    int i = 1;
    while (i + 1 < get_number_of_columns() && matrix::operator()(1, i + 1) < t) i++;
    double res = 0.0;
    double r = (*this)(t);
    res += (matrix::operator()(2, i + 1) + r) * (matrix::operator()(1, i + 1) - t) / 2;
    i++;

    while (i + 1 < get_number_of_columns() && matrix::operator()(1, i + 1) < T) {
        res += (matrix::operator()(2, i + 1) + matrix::operator()(2, i)) * (matrix::operator()(1, i + 1) - matrix::operator()(1, i)) / 2;
        i++;
    }
    r = (*this)(T);
    res += (r + matrix::operator()(2, i)) * (T - matrix::operator()(1, i)) / 2;

    return res;
}