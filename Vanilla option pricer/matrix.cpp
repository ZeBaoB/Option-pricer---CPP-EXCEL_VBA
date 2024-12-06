#include "pch.h"
#include "matrix.h"

// Assignment operator
matrix& matrix::operator=(const matrix& m) {
    if (this != &m) {
        M_ = m.M_;
        N_ = m.N_;
        data_ = m.data_;
    }
    return *this;
}

// Element access operators (1-based indexing)
  long double matrix::operator()(unsigned int l, unsigned int c) const {
    if (l < 1 || l > M_ || c < 1 || c > N_)
        throw std::out_of_range("Index out of bounds");
    return data_[l - 1][c - 1]; // Adjust for 1-based indexing
}

  long double& matrix::operator()(unsigned int l, unsigned int c) {
    if (l < 1 || l > M_ || c < 1 || c > N_)
        throw std::out_of_range("Index out of bounds");
    return data_[l - 1][c - 1]; // Adjust for 1-based indexing
}

// Fill a line with a given value (1-based indexing)
matrix& matrix::fill_line(unsigned int i,   long double alpha) {
    if (i < 1 || i > M_) throw std::out_of_range("Line index out of bounds");
    for (unsigned int j = 0; j < N_; ++j) data_[i - 1][j] = alpha;
    return *this;
}

// Fill a column with a given value (1-based indexing)
matrix& matrix::fill_column(unsigned int j,   long double alpha) {
    if (j < 1 || j > N_) throw std::out_of_range("Column index out of range");
    for (unsigned int i = 0; i < M_; ++i) data_[i][j - 1] = alpha;
    return *this;
}

// Add a scaled version of one row to another
matrix& matrix::add_line_to_line(unsigned int i1,   long double alpha, unsigned int i2) {
    if (i1 < 1 || i1 > M_ || i2 < 1 || i2 > M_)
        throw std::out_of_range("Line index out of bounds");
    for (unsigned int j = 0; j < N_; ++j) {
        data_[i1 - 1][j] += alpha * data_[i2 - 1][j];
    }
    return *this;
}

// Scalar multiplication (1-based indexing)
matrix operator*(const   long double k, const matrix& M1) {
    matrix m(M1.get_number_of_lines(), M1.get_number_of_columns());
    for (unsigned int i = 1; i <= m.get_number_of_lines(); ++i) {
        for (unsigned int j = 1; j <= m.get_number_of_columns(); ++j) {
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
    for (unsigned int i = 1; i <= m.get_number_of_lines(); ++i) {
        for (unsigned int j = 1; j <= m.get_number_of_columns(); ++j) {
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
    unsigned int nl = M1.get_number_of_lines();
    unsigned int nc = M2.get_number_of_columns();
    matrix m(nl, nc);

    for (unsigned int i = 1; i <= nl; ++i) {
        for (unsigned int j = 1; j <= nc; ++j) {
            m(i, j) = 0.0;
            for (unsigned int k = 1; k <= M1.get_number_of_columns(); ++k) {
                m(i, j) += M1(i, k) * M2(k, j);
            }
        }
    }
    return m;
}

// Overload << for matrix output (1-based indexing)
std::ostream& operator<<(std::ostream& st, const matrix& M) {
    unsigned int nl = M.get_number_of_lines();
    unsigned int nc = M.get_number_of_columns();
    st << "{";
    for (unsigned int i = 1; i <= nl; ++i) {
        for (unsigned int j = 1; j <= nc; ++j) {
            st << M(i, j) << (j == nc ? "" : ", ");
        }
        st << (i == nl ? "" : "\n ");
    }
    st << "}";
    return st;
}