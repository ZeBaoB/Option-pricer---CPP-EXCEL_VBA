#include "pch.h"
#include "matrix_plf.h"
#include <stdexcept>
#include <memory>


long double  matrix_plf::operator()(long double  t) const {
    unsigned int columns = get_number_of_columns();
    // The operator returns the first (last) value of returns if the   long double  t is lower than the lowest time (is higher than the highest time) in the tab
    if (t <= matrix::operator()(1, 1)) { return matrix::operator()(2, 1); }
    else if (t >= matrix::operator()(1, columns)) { return matrix::operator()(2, columns); }
    else {
        // Linear interpolation for t in [t_i; t_i+1]
        long double  t1, t2, r1, r2;
        for (unsigned int j = 1; j < columns; ++j) {
            t1 = matrix::operator()(1, j);
            t2 = matrix::operator()(1, j + 1);
            if (t >= t1 && t <= t2) {
                r1 = matrix::operator()(2, j);
                r2 = matrix::operator()(2, j + 1);
                r2 = matrix::operator()(2, j + 1);
                return r1 + (r2 - r1) * (t - t1) / (t2 - t1);
            }
        }
    }
    return matrix::operator()(2, columns);
}


matrix_plf matrix_plf::shift_value(long double  delta_r) {
    unsigned int numberOfPoints = get_number_of_columns();

    // Create 2 smart pointers for arrays of times and returns
    auto tab_t = std::make_unique< long double[]>(numberOfPoints);
    auto tab_r = std::make_unique< long double[]>(numberOfPoints);

    // Fill `tab_t` and `tab_r` with the current values
    for (unsigned int j = 1; j <= numberOfPoints; j++) {
        tab_t[j - 1] = matrix::operator()(1, j);          // Copy the first row (t values)
        tab_r[j - 1] = matrix::operator()(2, j) + delta_r; // Shift the second row (r values) by `delta_r`
    }

    // Create a new `matrix_plf` with the shifted values
    matrix_plf m(tab_t.get(), tab_r.get(), numberOfPoints);

    return m;
}


long double  matrix_plf::integral(long double  t, long double  T) const
{
    unsigned int i = 1;
    while (i + 1 < get_number_of_columns() && matrix::operator()(1, i + 1) < t) i++;
    long double  res = 0.0;
    long double  r = (*this)(t);
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