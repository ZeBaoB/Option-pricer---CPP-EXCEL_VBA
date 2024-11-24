#include "pch.h"
#include "matrix.h"
#include <iostream>

// Assignment operator
matrix& matrix::operator=(const matrix& m) {
    if (this == &m) {
        return *this;  // Check for self-assignment
    }
    delete[] data_;
    int nl = m.get_number_of_lines();
    int nc = m.get_number_of_columns();
    M_ = nl;
    N_ = nc;
    data_ = new double[nl * nc];
    std::copy(m.data_, m.data_ + M_ * N_, data_);
    return *this;
}

// Element access operators (1-based indexing)
double matrix::operator()(int l, int c) const {
    return data_[(l - 1) * N_ + (c - 1)];
}

double& matrix::operator()(int l, int c) {
    return data_[(l - 1) * N_ + (c - 1)];
}

// Fill a line with a given value (1-based indexing)
matrix& matrix::fill_line(int i, double alpha) {
    if (i < 1 || i > M_) {
        throw std::out_of_range("Line index out of range");
    }
    for (int j = 1; j <= N_; ++j) {
        (*this)(i, j) = alpha;
    }
    return *this;
}

// Fill a column with a given value (1-based indexing)
matrix& matrix::fill_column(int j, double alpha) {
    if (j < 1 || j > N_) {
        throw std::out_of_range("Column index out of range");
    }
    for (int i = 1; i <= M_; ++i) {
        (*this)(i, j) = alpha;
    }
    return *this;
}

// Add one line to another with a multiplier (1-based indexing)
matrix& matrix::add_line_to_line(int i1, double alpha, int i2) {
    if (i1 < 1 || i1 > M_ || i2 < 1 || i2 > M_) {
        throw std::out_of_range("Line index out of range");
    }
    for (int j = 1; j <= N_; ++j) {
        (*this)(i1, j) += alpha * (*this)(i2, j);
    }
    return *this;
}

matrix matrix::transpose()
{
    int nl = get_number_of_columns();
    int nc = get_number_of_lines();

    matrix m = matrix(nl, nc);
    for (int i = 1; i <= nl; ++i) { for (int j = 1; j <= nc; ++j) { m(i, j) = (*this)(j, i); } }

    return m;
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