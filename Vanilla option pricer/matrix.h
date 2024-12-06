#pragma once
#include "pch.h"
#include <iostream>
#include <vector>
/**
 * @file matrix.h
 * @brief Defines the `matrix` class, a general-purpose matrix class that serves as the parent for other specialized matrix types.
 *
 * @details
 * The `matrix` class provides basic functionality for matrix operations, including element access, modification,
 * and common operations like scaling, addition, and multiplication. It uses a 2D `std::vector` to store matrix data.
 *
 * @authors
 * - Bill Njonze
 * - Giorgio Enrico Maria Serva
 *
 * @date November 2024
 */

class matrix {
private:
	unsigned int M_; ///< Number of rows in the matrix.
	unsigned int N_; ///< Number of columns in the matrix.
	std::vector< std::vector<  long double > > data_; ///< Matrix data stored as a 2D vector.
public:
	/**
	 * @brief Constructs a matrix with specified dimensions.
	 *
	 * @param nl Number of rows.
	 * @param nc Number of columns.
	 *
	 * @details
	 * Initializes the matrix with all elements set to 0.0.
	 */
	matrix(unsigned int nl, unsigned int nc) : M_(nl), N_(nc), data_(nl, std::vector< long double>(nc, 0.0)) {}

	/**
	 * @brief Retrieves the number of rows in the matrix.
	 *
	 * @return The number of rows (`M_`).
	 */
	unsigned int get_number_of_lines() const { return M_; }

	/**
	* @brief Retrieves the number of columns in the matrix.
	*
	* @return The number of columns (`N_`).
	*/
	unsigned int get_number_of_columns() const { return N_; }

	/**
	 * @brief Overloads the assignment operator to copy another matrix.
	 *
	 * @param m The matrix to be copied.
	 * @return A reference to the updated matrix.
	 */
	matrix& operator=(const matrix& m);

	/**
	 * @brief Accesses an element of the matrix (read-only).
	 *
	 * @param l The row index (1-based).
	 * @param c The column index (1-based).
	 * @return The value of the matrix element at the specified position.
	 *
	 * @note Indices are assumed to start at 1 for this implementation.
	 */
	long double operator()(unsigned int l, unsigned int c) const;

	/**
	 * @brief Accesses and modifies an element of the matrix.
	 *
	 * @param l The row index (1-based).
	 * @param c The column index (1-based).
	 * @return A reference to the matrix element at the specified position.
	 *
	 * @note Indices are assumed to start at 1 for this implementation.
	 */
	long double& operator()(unsigned int l, unsigned int c);

	/**
	* @brief Fills a specified row with a given value.
	*
	* @param i The row index (1-based).
	* @param alpha The value to fill the row with.
	* @return A reference to the modified matrix.
	*/
	matrix& fill_line(unsigned int i, long double alpha);

	/**
	 * @brief Fills a specified column with a given value.
	 *
	 * @param j The column index (1-based).
	 * @param alpha The value to fill the column with.
	 * @return A reference to the modified matrix.
	 */
	matrix& fill_column(unsigned int j, long double alpha);

	/**
	 * @brief Adds a scaled row to another row.
	 *
	 * @param i1 The destination row index (1-based).
	 * @param alpha The scalar multiplier for the source row.
	 * @param i2 The source row index (1-based).
	 * @return A reference to the modified matrix.
	 *
	 * @details Performs the operation: `row[i1] += alpha * row[i2]`.
	 */
	matrix& add_line_to_line(unsigned int i1, long double alpha, unsigned int i2); // line i1 <-- line i1 + alpha * line i2
};

/**
 * @brief Scales a matrix by a scalar value.
 *
 * @param k The scalar multiplier.
 * @param M1 The matrix to be scaled.
 * @return A new matrix representing the result of the scaling operation.
 */
matrix operator*(const   long double k, const matrix& M1);

/**
 * @brief Adds two matrices element-wise.
 *
 * @param M1 The first matrix.
 * @param M2 The second matrix.
 * @return A new matrix representing the element-wise sum of `M1` and `M2`.
 */
matrix operator+(const matrix& M1, const matrix& M2);

/**
 * @brief Multiplies two matrices.
 *
 * @param M1 The first matrix.
 * @param M2 The second matrix.
 * @return A new matrix representing the result of the matrix multiplication.
 *
 * @details Assumes that the number of columns in `M1` matches the number of rows in `M2`.
 */
matrix operator*(const matrix& M1, const matrix& M2);

/**
 * @brief Outputs a matrix to an output stream.
 *
 * @param st The output stream (e.g., `std::cout`).
 * @param M The matrix to be displayed.
 * @return A reference to the output stream.
 */
std::ostream& operator<<(std::ostream& st, const matrix& M);