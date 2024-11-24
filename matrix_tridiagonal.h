#pragma once
#include "matrix.h"

class matrix_tridiagonal :
    public matrix
{
private:
	matrix::matrix;
public:

	matrix_tridiagonal(int dim) : matrix(dim, dim) {};
	matrix_tridiagonal(int dim, double* diag, double* diagSup, double* diagInf) : matrix(dim, dim) {
		matrix::operator()(1, 1) = diag[0];
		for (int i = 2; i <= dim; ++i) {
			matrix::operator()(i, i) = diag[i - 1];
			matrix::operator()(i - 1, i) = diagSup[i - 2];
			matrix::operator()(i, i - 1) = diagInf[i - 2];
		}
	}
	matrix inverse() const;
};

