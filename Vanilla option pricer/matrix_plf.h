#pragma once
#include "matrix.h"
//Matrix for partial linear function

class matrix_plf :
    public matrix
{
private:
	matrix::matrix;
public:
	// Constructor with 1-based indexing for input arrays `tab_t` and `tab_r`
	matrix_plf(int numberOfPoints) : matrix(2, numberOfPoints) {}
	matrix_plf(const double* tab_t, const double* tab_r, int numberOfPoints) : matrix(2, numberOfPoints) {
		// Directly access elements with 1-based indexing
		for (int j = 1; j <= numberOfPoints; ++j) {
			matrix::operator()(1, j) = tab_t[j - 1]; // Fill row 1 with `tab_t`
			matrix::operator()(2, j) = tab_r[j - 1]; // Fill row 2 with `tab_r`
		}
	}

	double operator()(double t) const;
	matrix_plf shift_value(double delta);
	double integral(double t, double T) const;
};