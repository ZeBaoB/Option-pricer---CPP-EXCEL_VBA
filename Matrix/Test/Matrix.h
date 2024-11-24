#include <iostream>
#include <stdexcept>
#include <memory> // Pour std::unique_ptr

class matrix {
private:
	int M_; // Number of lines
	int N_; // Number of columns
	std::unique_ptr<double[]> data_; // number stored in the matrix

	int index(int l, int c) const {  // Helper to convert indices
		return (l - 1) * N_ + (c - 1);
	}

public:
	// Constructeur
	matrix(int nl, int nc) : M_(nl), N_(nc), data_(std::make_unique<double[]>(nl* nc)) {
		std::fill(data_.get(), data_.get() + nl * nc, 0.0); // Initialiser à 0
	}


	int get_number_of_lines() const { return M_; }
	int get_number_of_columns() const { return N_; }
	matrix& operator=(const matrix& other);

	// element access operator for reading
	double operator()(int l, int c) const;
	// element access operator for modification
	double& operator()(int l, int c);

	matrix& fill_line(int i, double alpha);
	matrix& fill_column(int j, double alpha);
	matrix& add_line_to_line(int i1, double alpha, int i2); // line i1 <-- line i1 + alpha * line i2

};

matrix operator*(const double k, const matrix& M1);
matrix operator+(const matrix& M1, const matrix& M2);
matrix operator*(const matrix& M1, const matrix& M2);
std::ostream& operator<<(std::ostream& st, const matrix& M);


class matrix_plf : public matrix {
private:
	matrix::matrix;

public:
	// Constructor with 1-based indexing for input arrays `tab_t` and `tab_r`
	matrix_plf(int numberOfPoints) : matrix(2, numberOfPoints) {}
	matrix_plf(double* tab_t, double* tab_r, int numberOfPoints) : matrix(2, numberOfPoints) {
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

class matrix_tridiagonal : public matrix {
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

class mesh_matrix : public matrix {
private:
	matrix::matrix;
	double T_;
	double s0_;
	double K_;
	double sinf_;
public:
	mesh_matrix(double T,
		double s0,
		double K,
		int spot_mest_parameter,
		int time_mesh_parameter,
		double sinf) : matrix(time_mesh_parameter, spot_mest_parameter) {
		T_ = T;
		s0_ = s0;
		K_ = K;
		sinf_ = sinf;
	}

	double get_T() { return T_; }
	double get_s0() { return s0_; }
	double get_K() { return K_; }
	double get_sinf() { return sinf_; }		

	void solve(bool call, bool european, matrix_plf r_table, double sigma);
	double retrieve_OptionValue(double s0, int line = 1);
	double retrieve_delta(double s0);
	double retrieve_gamma(double s0);
	double retrieve_theta(double s0); // % t

	double retrieve_rho(bool call, bool european, matrix_plf r_table, double sigma, double s0); // %r
	double retrieve_vega(bool call, bool european, matrix_plf r_table, double sigma, double s0); // %sigma

};
