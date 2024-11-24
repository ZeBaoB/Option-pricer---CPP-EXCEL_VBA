#pragma once
#include "matrix.h"
#include "matrix_plf.h"
#include "matrix_tridiagonal.h"
#include <cmath>

class mesh_matrix :
	public matrix
{
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