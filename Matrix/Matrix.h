#pragma once
#include <iostream>
#include <vector>

class matrix {
private:
	int M_; // Number of lines
	int N_; // Number of columns
	std::vector< std::vector<long double> > data_; // number stored in the matrix
public:
	// Constructor
	matrix(int nl, int nc) : M_(nl), N_(nc), data_(nl, std::vector<long double>(nc, 0.0)) {}

	int get_number_of_lines() const { return M_; }
	int get_number_of_columns() const { return N_; }
	matrix& operator=(const matrix& m);

	// element access operator for reading
	long double operator()(int l, int c) const;
	// element access operator for modification
	long double& operator()(int l, int c);

	matrix& fill_line(int i, long double alpha);
	matrix& fill_column(int j, long double alpha);
	matrix& add_line_to_line(int i1, long double alpha, int i2); // line i1 <-- line i1 + alpha * line i2
};

matrix operator*(const long double k, const matrix& M1);
matrix operator+(const matrix& M1, const matrix& M2);
matrix operator*(const matrix& M1, const matrix& M2);
std::ostream& operator<<(std::ostream& st, const matrix& M);