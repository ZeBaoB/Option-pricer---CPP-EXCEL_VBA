#pragma once
#include "pch.h"
#include "matrix.h"

/**
 * @file matrix_tridiagonal.h
 * @brief Defines the `matrix_tridiagonal` class, a specialized class for tridiagonal matrices.
 *
 * @details
 * The `matrix_tridiagonal` class extends the base `matrix` class and provides functionality
 * for constructing and manipulating tridiagonal matrices. These matrices have non-zero elements
 * only on the main diagonal, the diagonal above it, and the diagonal below it.
 *
 * @authors
 * - Bill Njonze
 * - Giorgio Enrico Maria Serva
 *
 * @date November 2024
 */
class matrix_tridiagonal :
	public matrix
{
private:
	matrix::matrix;
public:

	/**
	 * @brief Constructs a tridiagonal matrix with specified dimensions.
	 *
	 * @param dim The dimension of the matrix (number of rows and columns).
	 *
	 * @details
	 * Initializes the matrix as a square matrix of size `dim x dim` with all elements set to zero.
	 */
	matrix_tridiagonal(int dim) : matrix(dim, dim) {};

	/**
	 * @brief Constructs a tridiagonal matrix with specified diagonals.
	 *
	 * @param dim The dimension of the matrix (number of rows and columns).
	 * @param diag Pointer to an array containing the main diagonal elements (length `dim`).
	 * @param diagSup Pointer to an array containing the elements of the superdiagonal (length `dim - 1`).
	 * @param diagInf Pointer to an array containing the elements of the subdiagonal (length `dim - 1`).
	 *
	 * @details
	 * Populates the matrix such that:
	 * - The main diagonal is filled with elements from `diag`.
	 * - The superdiagonal (above the main diagonal) is filled with elements from `diagSup`.
	 * - The subdiagonal (below the main diagonal) is filled with elements from `diagInf`.
	 */
	matrix_tridiagonal(int dim, double* diag, double* diagSup, double* diagInf) : matrix(dim, dim) {
		matrix::operator()(1, 1) = diag[0];
		for (int i = 2; i <= dim; ++i) {
			matrix::operator()(i, i) = diag[i - 1];
			matrix::operator()(i - 1, i) = diagSup[i - 2];
			matrix::operator()(i, i - 1) = diagInf[i - 2];
		}
	}

	/**
	 * @brief Solves the linear system using the Gauss-Seidel method.
	 *
	 * @param b The matrix representing the right-hand side of the system.
	 * @param x0 The initial guess for the solution.
	 * @return A matrix representing the solution to the linear system.
	 *
	 * @details
	 * This method applies the Gauss-Seidel iterative algorithm to solve the linear system
	 * represented by the tridiagonal matrix and the vector `b`.
	 */
	matrix gauss_saidel(matrix b, matrix x0) const;
	
	/**
	 * @brief Computes the inverse of the tridiagonal matrix.
	 *
	 * @return A matrix representing the inverse of the tridiagonal matrix.
	 *
	 * @details
	 * This method calculates the inverse of the matrix if it is invertible.
	 */
	matrix inverse() const;

};