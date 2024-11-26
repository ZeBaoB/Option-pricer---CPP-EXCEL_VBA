#pragma once
#include "pch.h"
#include "matrix.h"

/**
 * @file matrix_plf.h
 * @brief Defines the `matrix_plf` class, which represents a matrix for partial linear functions.
 *
 * The class extends the base `matrix` class and provides specific functionality for handling
 * partial linear functions, such as interpolation, value shifting, and integral computation.
 *
 * @details
 * - Each row of the matrix corresponds to a different data series:
 *   - Row 1: Time points (`tab_t`).
 *   - Row 2: Corresponding values (`tab_r`).
 * - Assumes 1-based indexing for input arrays.
 *
 * @authors
 * - Bill Njonze
 * - Giorgio Enrico Maria Serva
 *
 * @date November 2024
 */
class matrix_plf :
    public matrix
{
private:
	matrix::matrix;
public:
	/**
	 * @brief Constructs a `matrix_plf` object with a specified number of points.
	 *
	 * @param numberOfPoints The number of points in the matrix.
	 * Initializes a 2-row matrix with the given number of columns.
	 */
	matrix_plf(int numberOfPoints) : matrix(2, numberOfPoints) {}

	/**
	 * @brief Constructs a `matrix_plf` object from input arrays with time-value pairs.
	 *
	 * @param tab_t Pointer to the array of time points (1-based indexing assumed).
	 * @param tab_r Pointer to the array of corresponding values (1-based indexing assumed).
	 * @param numberOfPoints The number of points in the input arrays.
	 *
	 * @details
	 * Populates the matrix such that:
	 * - Row 1 contains elements from `tab_t`.
	 * - Row 2 contains elements from `tab_r`.
	 */
	matrix_plf(const double* tab_t, const double* tab_r, int numberOfPoints) : matrix(2, numberOfPoints) {
		// Directly access elements with 1-based indexing
		for (int j = 1; j <= numberOfPoints; ++j) {
			matrix::operator()(1, j) = tab_t[j - 1]; // Fill row 1 with `tab_t`
			matrix::operator()(2, j) = tab_r[j - 1]; // Fill row 2 with `tab_r`
		}
	}

	/**
	 * @brief Computes the value at a specific time `t` using interpolation.
	 *
	 * @param t The time point for which to compute the value.
	 * @return The interpolated value if `t` is between two dates in the table.
	 * Returns:
	 * - The value at the first time point if `t` is before the first date.
	 * - The value at the last time point if `t` is after the last date.
	 */
	double operator()(double t) const;

	/**
	 * @brief Shifts all values in the matrix by a specified delta.
	 *
	 * @param delta The value by which to shift the risk-free rate.
	 * @return A new `matrix_plf` object with all values shifted by `delta`.
	 */
	matrix_plf shift_value(double delta);
	
	/**
	 * @brief Computes the discount rate between two time points.
	 *
	 * @param t The start time.
	 * @param T The end time.
	 * @return The computed discount rate between `t` and `T`.
	 */
	double integral(double t, double T) const;
};