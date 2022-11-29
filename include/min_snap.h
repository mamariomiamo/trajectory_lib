#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iosqp.hpp"
#include <iostream>

namespace MIN_SNAP
{
    struct poly_traj
    {
        // coeff has 3 entries for 3 dimensions
        // each entry has k_segment of columns with each column as the the polynomial coefficients [p_0, p_1, ..., p_n]'
        // p(t) = p_0 + p_1 * t + p_2 * t ^ 2 + ... + p_n * t ^ n
        std::vector<Eigen::MatrixXd> coeff;
        std::vector<double> time_duration;

        // used just to store values for single axis
        Eigen::VectorXd x_coeff;
        Eigen::VectorXd y_coeff;
        Eigen::VectorXd z_coeff;
        Eigen::MatrixXd x_coeff_mat;
        Eigen::MatrixXd y_coeff_mat;
        Eigen::MatrixXd z_coeff_mat;
    };

    class min_snap_traj
    {
    public:
        min_snap_traj(int m_dimension, int k_segment, int n_order, int r_derivative_order, bool close_form);

        min_snap_traj(const std::vector<Eigen::Vector3d> &waypt_list, int m_dimension, bool close_form);

        min_snap_traj(const std::vector<Eigen::Vector3d> &waypt_list, int m_dimension, int n_order, int r_derivative_order, bool close_form);

        ~min_snap_traj() = default;

        // helper functions

        int factorial(int x);

        double computeDistance(const std::vector<Eigen::Vector3d> &waypt_list, std::vector<double> &segment_length);

        std::vector<double> computeT(const std::vector<Eigen::Vector3d> &waypt_list, const Eigen::Vector3d &vel_max);

        std::vector<double> computeT(const std::vector<Eigen::Vector3d> &waypt_list, double total_time);

        void evaluatePoly(Eigen::VectorXd &row_vector, int derivative_order, double time, int polynomial_order);

        void computeHessian(Eigen::MatrixXd &hessian, int hessian_dim, int n_poly_order, int r_derivate_order, int k_segment, const std::vector<double> &time_vector);

        void computeEquality(Eigen::MatrixXd &Aeq, int Aeq_row, int Aeq_col, Eigen::VectorXd &beq, int beq_dim, int n_poly_order, int k_segment, const std::vector<double> &waypoint_list, const std::vector<double> &time_vector);

        void computeAeq(Eigen::MatrixXd &Aeq, int n_poly_order, int k_segment, const std::vector<double> &time_vector);

        void computeBeq(Eigen::VectorXd &beq, int k_segment, const std::vector<double> &waypoint_list);

        void computeDerivateMappingMatrix(Eigen::MatrixXd &mapping_A, int k_segment, const std::vector<double> &time_vector);

        Eigen::VectorXd optimize(const Eigen::MatrixXd &hessian_Q, const Eigen::MatrixXd &equality_A, Eigen::VectorXd &equality_b);

        // polynomial coefficients will be computed using constrained QP and stored in coeff
        void computeCoeffConstrainedQP(Eigen::VectorXd &coeff, const Eigen::MatrixXd &hessian_Q, int n_poly_order, int k_segment, const std::vector<double> &waypt_single_axis, const std::vector<double> &time_vector);

        // polynomial coefficients will be computed using closed-form method and stored in coeff
        void computeCoeffCloseFormQPSingleAxis(Eigen::VectorXd &coeff, const Eigen::MatrixXd &hessian_Q, int n_poly_order, int k_segment, const std::vector<double> &waypt_single_axis, const std::vector<double> &time_vector);

        void computeCoeffCloseFormQPThreeAxis(std::vector<Eigen::MatrixXd> &coeff, const Eigen::MatrixXd &hessian_Q, int n_poly_order, int k_segment, const std::vector<Eigen::Vector3d> &waypt_list, const std::vector<double> &time_vector);

        void coeffMatReshape(const Eigen::MatrixXd &coeff_all, std::vector<Eigen::MatrixXd> &coeff);

        poly_traj computeTrajectory(const std::vector<Eigen::Vector3d> &waypt_list);

    private:
        int m_dimension, k_segment;
        int n_order = 5;
        int r_derivative_order = 4;
        int v_norm = 2;
        bool close_form;
        std::vector<Eigen::Vector3d> waypt_list;
        std::vector<double> time_vector;
        poly_traj trajectory;
    };

} // namespace MIN_SNAP_QP