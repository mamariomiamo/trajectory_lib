#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "iosqp.hpp"

using namespace Eigen;
using namespace std;

namespace MIN_SNAP_QP
{
    class min_snap_traj
    {
    public:
        min_snap_traj();
        ~min_snap_traj() = default;

        // helper functions

        int factorial(int x);

        vector<double> computeT(const vector<Vector3d> &waypt_list, const Vector3d &vel_max);

        void evaluatePoly(VectorXd &row_vector, int derivative_order, double time, int polynomial_order);

        void computeHessian(MatrixXd &hessian, int hessian_dim, int n_poly_order, int r_derivate_order, int k_segment, const vector<double> &time_vector);

        void computeEquality(MatrixXd &Aeq, int Aeq_row, int Aeq_col, VectorXd &beq, int beq_dim, int n_poly_order, int k_segement, const vector<double> &waypoint_list, const vector<double> &time_vector);

        void computeAeq(MatrixXd &Aeq, int n_poly_order, int k_segement, const vector<double> &time_vector);

        void computeBeq(VectorXd &beq, int k_segement, const vector<double> &waypoint_list);

        VectorXd optimize(const MatrixXd &hessian_Q, const MatrixXd &equality_A, VectorXd &equality_b);

        // polynomial coefficients will be computed using constrained QP and stored in coeff        
        void computeCoeffConstrainedQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector);

        // polynomial coefficients will be computed using closed-form method and stored in coeff
        void computeCoeffCloseFormQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector);

    private:
        int m_dimension, k_segment, n_order, r_derivative_order;
        bool close_form;
    };

} // namespace MIN_SNAP_QP