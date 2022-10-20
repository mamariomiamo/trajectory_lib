#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iosqp.hpp"
#include <iostream>

#define log(x) std::cout << x << std::endl;

using namespace Eigen;
using namespace std;

namespace MIN_SNAP
{
    struct poly_traj
    {
        // coeff has 3 entries for 3 dimensions
        // each entry has k_segment of columns with each column as the the polynomial coefficients [p_0, p_1, ..., p_n]'
        // p(t) = p_0 + p_1 * t + p_2 * t ^ 2 + ... + p_n * t ^ n
        vector<MatrixXd> coeff;
        vector<double> time_duration;

        // used just to store values for single axis
        VectorXd x_coeff;
        VectorXd y_coeff;
        VectorXd z_coeff;
        MatrixXd x_coeff_mat;
        MatrixXd y_coeff_mat;
        MatrixXd z_coeff_mat;
    };

    class min_snap_traj
    {
    public:
        min_snap_traj(int m_dimension, int k_segment, int n_order, int r_derivative_order, bool close_form);

        min_snap_traj(const vector<Vector3d> &waypt_list, int m_dimension, bool close_form);

        min_snap_traj(const vector<Vector3d> &waypt_list, int m_dimension, int n_order, int r_derivative_order, bool close_form);

        ~min_snap_traj() = default;

        // helper functions

        int factorial(int x);

        double computeDistance(const vector<Vector3d> &waypt_list, vector<double> &segment_length);

        vector<double> computeT(const vector<Vector3d> &waypt_list, const Vector3d &vel_max);

        vector<double> computeT(const vector<Vector3d> &waypt_list, double total_time);

        void evaluatePoly(VectorXd &row_vector, int derivative_order, double time, int polynomial_order);

        void computeHessian(MatrixXd &hessian, int hessian_dim, int n_poly_order, int r_derivate_order, int k_segment, const vector<double> &time_vector);

        void computeEquality(MatrixXd &Aeq, int Aeq_row, int Aeq_col, VectorXd &beq, int beq_dim, int n_poly_order, int k_segment, const vector<double> &waypoint_list, const vector<double> &time_vector);

        void computeAeq(MatrixXd &Aeq, int n_poly_order, int k_segment, const vector<double> &time_vector);

        void computeBeq(VectorXd &beq, int k_segment, const vector<double> &waypoint_list);

        void computeDerivateMappingMatrix(MatrixXd &mapping_A, int k_segment, const vector<double> &time_vector);

        VectorXd optimize(const MatrixXd &hessian_Q, const MatrixXd &equality_A, VectorXd &equality_b);

        // polynomial coefficients will be computed using constrained QP and stored in coeff
        void computeCoeffConstrainedQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector);

        // polynomial coefficients will be computed using closed-form method and stored in coeff
        void computeCoeffCloseFormQPSingleAxis(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector);

        void computeCoeffCloseFormQPThreeAxis(vector<MatrixXd> &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<Vector3d> &waypt_list, const vector<double> &time_vector);

        void coeffMatReshape(const MatrixXd &coeff_all, vector<MatrixXd> &coeff);

        poly_traj computeTrajectory(const vector<Vector3d> &waypt_list);

    private:
        int m_dimension, k_segment;
        int n_order = 5;
        int r_derivative_order = 4;
        int v_norm = 2;
        bool close_form;
        vector<Vector3d> waypt_list;
        vector<double> time_vector;
        poly_traj trajectory;
    };

} // namespace MIN_SNAP_QP