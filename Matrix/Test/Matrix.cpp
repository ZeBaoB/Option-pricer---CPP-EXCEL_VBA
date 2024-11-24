#include "Matrix.h"
#include <stdexcept>
#include <cmath>

// Assignment operator
matrix& matrix::operator=(const matrix& other) {
    if (this == &other) {
        return *this;  // Check for self-assignment
    }
    M_ = other.M_;
    N_ = other.N_;
    data_ = std::make_unique<double[]>(M_ * N_);
    std::copy(other.data_.get(), other.data_.get() + M_ * N_, data_.get());
}

// Element access operators (1-based indexing)
double matrix::operator()(int l, int c) const {
    if (l < 1 || l > M_ || c < 1 || c > N_) throw std::out_of_range("Matrix index out of bounds");
    return data_[index(l, c)];
}

double& matrix::operator()(int l, int c) {
    if (l < 1 || l > M_ || c < 1 || c > N_) throw std::out_of_range("Matrix index out of bounds");
    return data_[index(l, c)];
}

// Fill a line with a given value (1-based indexing)
matrix& matrix::fill_line(int i, double alpha) {
    if (i < 1 || i > M_) throw std::out_of_range("Line index out of bounds");
    for (int j = 1; j <= N_; ++j) (*this)(i, j) = alpha;
    return *this;
}

// Fill a column with a given value (1-based indexing)
matrix& matrix::fill_column(int j, double alpha) {
    if (j < 1 || j > N_) throw std::out_of_range("Column index out of range");
    for (int i = 1; i <= M_; ++i) (*this)(i, j) = alpha;
    return *this;
}

// Add one line to another with a multiplier (1-based indexing)
matrix& matrix::add_line_to_line(int i1, double alpha, int i2) {
    if (i1 < 1 || i1 > M_ || i2 < 1 || i2 > M_) throw std::out_of_range("Line index out of range");
    for (int j = 1; j <= N_; ++j) (*this)(i1, j) += alpha * (*this)(i2, j);
    return *this;
}

// Scalar multiplication (1-based indexing)
matrix operator*(const double k, const matrix& M1) {
    matrix m(M1.get_number_of_lines(), M1.get_number_of_columns());
    for (int i = 1; i <= m.get_number_of_lines(); ++i) {
        for (int j = 1; j <= m.get_number_of_columns(); ++j) {
            m(i, j) = k * M1(i, j);
        }
    }
    return m;
}

// Matrix addition (1-based indexing)
matrix operator+(const matrix& M1, const matrix& M2) {
    if (M1.get_number_of_lines() != M2.get_number_of_lines() || M1.get_number_of_columns() != M2.get_number_of_columns()) {
        throw std::invalid_argument("Dimension mismatch for matrix addition");
    }

    matrix m(M1.get_number_of_lines(), M1.get_number_of_columns());
    for (int i = 1; i <= m.get_number_of_lines(); ++i) {
        for (int j = 1; j <= m.get_number_of_columns(); ++j) {
            m(i, j) = M1(i, j) + M2(i, j);
        }
    }
    return m;
}

// Matrix multiplication (1-based indexing)
matrix operator*(const matrix& M1, const matrix& M2) {
    if (M1.get_number_of_columns() != M2.get_number_of_lines()) {
        throw std::invalid_argument("Dimension mismatch for matrix multiplication");
    }
    int nl = M1.get_number_of_lines();
    int nc = M2.get_number_of_columns();
    matrix m(nl, nc);

    for (int i = 1; i <= nl; ++i) {
        for (int j = 1; j <= nc; ++j) {
            m(i, j) = 0.0;
            for (int k = 1; k <= M1.get_number_of_columns(); ++k) {
                m(i, j) += M1(i, k) * M2(k, j);
            }
        }
    }
    return m;
}

// Overload << for matrix output (1-based indexing)
std::ostream& operator<<(std::ostream& st, const matrix& M) {
    int nl = M.get_number_of_lines();
    int nc = M.get_number_of_columns();
    st << "{";
    for (int i = 1; i <= nl; ++i) {
        for (int j = 1; j <= nc; ++j) {
            st << M(i, j) << (j == nc ? "" : ", ");
        }
        st << (i == nl ? "" : "\n ");
    }
    st << "}";
    return st;
}

double matrix_plf::operator()(double t) const {
    int columns = get_number_of_columns();
    // The operator returns the first (last) value of returns if the double t is lower than the lowest time (is higher than the highest time) in the tab
    if (t <= matrix::operator()(1, 1)) { return matrix::operator()(2, 1); }
    else if (t >= matrix::operator()(1, columns)) { return matrix::operator()(2, columns); }
    else{
        // Linear interpolation for t in [t_i; t_i+1]
        double t1, t2, r1, r2;
        for (int j = 1; j < columns; ++j) {
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
    while (i + 1 < get_number_of_columns() && matrix::operator()(1, i+1) < t) i++;
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

matrix matrix_tridiagonal::inverse() const {
    int dim = get_number_of_lines();

    // Copy of the current matrix to work on for inversion
    matrix m = matrix(dim, dim);
    for (int i = 1; i <= dim; ++i) { for (int j = 1; j <= dim; ++j) { m(i, j) = (*this)(i, j); } }

    // Initialize the identity matrix to store the inverse
    matrix m_1(dim, dim);
    for (int i = 1; i <= dim; ++i) {
        m_1(i, i) = 1.0;
    }

    // Lower triangular elimination
    for (int i = 1; i < dim; ++i) {
        double mii = m(i, i);
        double alpha = m(i + 1, i);

        m.add_line_to_line(i + 1, -alpha / mii, i);
        m_1.add_line_to_line(i + 1, -alpha / mii, i);
    }

    // Upper triangular elimination
    for (int i = dim; i > 1; --i) {
        double mii = m(i, i);
        double alpha = m(i - 1, i);

        m.add_line_to_line(i - 1, -alpha / mii, i);
        m_1.add_line_to_line(i - 1, -alpha / mii, i);
    }

    // Normalize diagonal to 1
    for (int i = 1; i <= dim; ++i) {
        double mii = m(i, i);
        for (int j = 1; j <= dim; ++j) {
            m_1(i, j) /= mii;
        }
    }
    return m_1;
};

void mesh_matrix::solve(bool call, bool european, matrix_plf r_table, double sigma)
{
    int dim = get_number_of_columns();
    int nl = get_number_of_lines();
    matrix Vi(dim, 1);
    matrix Vi_1(dim, 1);
    double r;
    double ak;
    double bk;
    double ck;
    double sk;
    double dt = T_ / (nl - 1);
    double S_K; // S-K or K - S

    if(call) {
        fill_column(1, 0.0);
        for (int j = 2; j <= dim; j++) {
            S_K = get_sinf() * (j - 1) / (dim - 1) - get_K();
            matrix::operator()(nl, j) = (S_K < 0.0) ? 0.0 : S_K;
        }
        for (int i = 1; i <= nl; i++) {
            matrix::operator()(i, dim) = get_sinf() - get_K() * std::exp(-r_table.integral( dt * (i - 1), get_T()));
        }
    }
    else {
        fill_column(dim, 0.0);
        for (int j = 2; j <= dim; j++) {
            S_K = get_K() - get_sinf() * (j - 1) / (dim - 1);
            matrix::operator()(nl, j) = (S_K < 0.0) ? 0.0 : S_K;
        }
        for (int i = 1; i <= nl; i++) {
            matrix::operator()(i, 1) = get_K() * std::exp(-r_table.integral( dt * (i - 1), get_T()));
        }
    }
    

    for (int i = nl; i >=2; --i) {
        for (int k = 1; k <= dim; ++k) {
            Vi(k, 1) = matrix::operator()(i, k);
        }

        r = (r_table( T_ * (i-1) / nl ) + r_table(T_ * i / nl))/2;

        matrix_tridiagonal Lmatrix(dim);
        matrix_tridiagonal Rmatrix(dim);
        /*
        Lmatrix(1, 1) = (matrix::operator()(i - 1, 1) > 0) ? 1 / matrix::operator()(i - 1, 1) : 1;
        Rmatrix(1, 1) = (Vi(1, 1) > 0.0) ? 1.0 / Vi(1, 1) : 1;
        Lmatrix(dim, dim) = (matrix::operator()(i - 1, dim) > 0) ? 1 / matrix::operator()(i - 1, dim) : 1;
        Rmatrix(dim, dim) = (Vi(dim, 1) > 0.0) ? 1.0 / Vi(dim, 1) : 1;
        */
        Lmatrix(1, 1) = 1.0;
        Rmatrix(1, 1) = 1.0;
        Lmatrix(dim, dim) = 1.0;
        Rmatrix(dim, dim) = 1.0;
        
        for (int k = 2; k <= dim - 1; k++) {
            sk = 1.0 * (k-1);
            ak = dt / 4 * (sigma * sigma* sk * sk - r * sk);
            bk = - dt / 2 * (sigma * sigma * sk * sk + r);
            ck = dt / 4 * (sigma * sigma * sk * sk + r * sk);

            Lmatrix(k, k - 1) = -ak;
            Lmatrix(k, k) = 1 - bk;
            Lmatrix(k, k + 1) = -ck;

            Rmatrix(k, k - 1) = ak;
            Rmatrix(k, k) = 1 + bk;
            Rmatrix(k, k + 1) = ck;
            
        }

        Vi_1 = Lmatrix.inverse() * Rmatrix * Vi;
        
        if (european){
            for (int k = 2; k <= dim-1; ++k) {
                matrix::operator()( i-1, k) = Vi_1( k, 1);
            }
        }
        else if (call) {
            for (int k = 2; k <= dim - 1; ++k) {
                matrix::operator()(i - 1, k) = (Vi_1(k, 1)>(get_sinf()*(k-1)/(dim - 1) - get_K())) ? Vi_1(k, 1) : (get_sinf() * (k - 1) / (dim - 1) - get_K());
            }
        }
        else {
            for (int k = 2; k <= dim - 1; ++k) {
                matrix::operator()(i - 1, k) = (Vi_1(k, 1) > (get_K() - get_sinf() * (k - 1) / (dim - 1))) ? Vi_1(k, 1) : (get_K() - get_sinf() * (k - 1) / (dim - 1));
            }
        }
    }
    //matrix::operator()(1, dim) = get_sinf() - get_K() * std::exp(-r_table.integral(0.0, get_T()));
}

double mesh_matrix::retrieve_OptionValue(double s0, int line)
{
    double j;
    int j1;
    int j2;
    j = 1.0 + s0 * (get_number_of_columns() - 1) / get_sinf();
    j1 = j;
    j2 = j + 1;

    return matrix::operator()(line, j1) + (j-j1) * (matrix::operator()(line, j2) - matrix::operator()(line, j1));
}

double mesh_matrix::retrieve_delta(double s0)
{
    double ds = get_sinf()/(get_number_of_columns() - 1);
    return (retrieve_OptionValue(s0 + ds) - retrieve_OptionValue(s0 - ds) ) / 2 / ds;
}

double mesh_matrix::retrieve_gamma(double s0)
{
    double ds = get_sinf() / (get_number_of_columns() - 1);
    return (retrieve_OptionValue(s0 + ds) + retrieve_OptionValue(s0 - ds) - 2 * retrieve_OptionValue(s0)) / ds / ds;
}

double mesh_matrix::retrieve_theta(double s0)
{
    return (retrieve_OptionValue(s0, 2) - retrieve_OptionValue(s0)) / get_T() * (get_number_of_lines() - 1);
}

double mesh_matrix::retrieve_rho(bool call, bool european, matrix_plf r_table, double sigma, double s0)
{
    double price1 = retrieve_OptionValue(s0);
    matrix_plf new_rtable = r_table.shift_value(0.001);
    solve(call, european, new_rtable, sigma);
    double price2 = retrieve_OptionValue(s0);

    solve(call, european, r_table, sigma);

    return (price2-price1)/0.001;
}

double mesh_matrix::retrieve_vega(bool call, bool european, matrix_plf r_table, double sigma, double s0)
{
    double price1 = retrieve_OptionValue(s0);
    solve(call, european, r_table, sigma + 0.001);
    double price2 = retrieve_OptionValue(s0);

    solve(call, european, r_table, sigma);

    return (price2 - price1) / 0.001;
}
