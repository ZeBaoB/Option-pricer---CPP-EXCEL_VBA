#pragma once
#include <iostream>

class matrix {
private:
	int M_; // Number of lines
	int N_; // Number of columns
	double* data_; // number stored in the matrix
public:
	matrix(int nl, int nc) {
		M_ = nl;
		N_ = nc;
		data_ = new double[nl * nc] {0};  // Dynamically allocate memory
	}

	int get_number_of_lines() const { return M_; }
	int get_number_of_columns() const { return N_; }
	matrix& operator=(const matrix& m);

	// element access operator for reading
	double operator()(int l, int c) const;
	// element access operator for modification
	double& operator()(int l, int c);

	matrix& fill_line(int i, double alpha);
	matrix& fill_column(int j, double alpha);
	matrix& add_line_to_line(int i1, double alpha, int i2); // line i1 <-- line i1 + alpha * line i2

	matrix transpose();
};

matrix operator*(const double k, const matrix& M1);
matrix operator+(const matrix& M1, const matrix& M2);
matrix operator*(const matrix& M1, const matrix& M2);
std::ostream& operator<<(std::ostream& st, const matrix& M);

