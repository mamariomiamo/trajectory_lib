#include "matplotlibcpp.h"
#include "min_snap.h"
#include "iosqp.hpp"
#include <iostream>
// #include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/Sparse>
#include <vector>
#include <iterator>

#define log(x) std::cout << x << std::endl;

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

int factorial(int x)
{
    int result = 1;
    for (int i = 1; i <= x; i++)
    {
        result *= i;
    }

    return result;
}

// example use of osqp to solve the quadratic programming problem as matlab example
// https://www.mathworks.com/help/optim/ug/quadprog.html#d124e140499
// f(x) = 1/2 * x^T H x + f^T x
// s.t. lb <= x <= ub
void testQuadProg()
{
    const double kInfinity = std::numeric_limits<double>::infinity();
    SparseMatrix<double> objective_matrix(2, 2);
    const Triplet<double> kTripletsP[] = {
        {0, 0, 1}, {1, 0, -1}, {0, 1, -1}, {1, 1, 2}};
    objective_matrix.setFromTriplets(std::begin(kTripletsP),
                                     std::end(kTripletsP));

    SparseMatrix<double> constraint_matrix(3, 2);
    const Triplet<double> kTripletsA[] = {
        {0, 0, 1}, {1, 0, -1}, {2, 0, 2}, {0, 1, 1}, {1, 1, 2}, {2, 1, 1}};
    constraint_matrix.setFromTriplets(std::begin(kTripletsA),
                                      std::end(kTripletsA));

    log("constraint is ");
    log(constraint_matrix);
    VectorXd objective_vector(2);
    objective_vector << -2, -6;

    VectorXd upper_bound(3);
    upper_bound << 2.0, 2.0, 3.0;

    VectorXd lower_bound(3);
    lower_bound << -kInfinity, -kInfinity, -kInfinity;

    IOSQP solver;
    c_int flag = solver.setMats(objective_matrix, objective_vector, constraint_matrix, lower_bound, upper_bound, 1e-3, 1e-3);
    if (flag != 0)
    {
        log("problem non-convex");
    }
    solver.solve();
    // c_int status = solver.getStatus();
    //  std::cout << "STATUS: " << status << std::endl;
    // q_ << 1.0, 1.0;

    log(solver.getPrimalSol());
}

double computeDistance(const vector<Vector3d> &waypt_list, vector<double> &segment_length);

vector<double> computeT(const vector<Vector3d> &waypt_list, const Vector3d &vel_max);

vector<double> computeT(const vector<Vector3d> &waypt_list, double total_time);

void computeHessian(MatrixXd &hessian, int hessian_dim, int n_poly_order, int r_derivate_order, int k_segment, const vector<double> &time_vector);

void computeEquality(MatrixXd &Aeq, int Aeq_row, int Aeq_col, VectorXd &beq, int beq_dim, int n_poly_order, int k_segment, const vector<double> &waypoint_list, const vector<double> &time_vector);

void computeAeq(MatrixXd &Aeq, int n_poly_order, int k_segment, const vector<double> &time_vector);

void computeBeq(VectorXd &beq, int k_segment, const vector<double> &waypoint_list);

/**
 *   \brief evaluatePoly.
 *
 *   Return a row vector according to the derivative order to evaluate, polynomial order and timestamp
 *
 *   \param derivative_order
 *   \param time
 *   \param polynomial_order
 *   \return a row vector
 *     e.g. [1 t t^2 t^3 ... t^n] for position, derivative_order = 0
 *          [0 1 2*t 3*t^2 ... n*t^(n-1)] for velocity, derivative_order = 1
 *
 **/
void evaluatePoly(VectorXd &row_vector, int derivative_order, double time, int polynomial_order);

VectorXd optimize(const MatrixXd &hessian_Q, const MatrixXd &equality_A, VectorXd &equality_b);
// void optimize(const MatrixXd &hessian_Q, const MatrixXd &equality_A, VectorXd &equality_b);

void computeDerivateMappingMatrix(MatrixXd &mapping_A, int k_segment, const vector<double> &time_vector);

void computeCoeffConstrainedQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector);

void computeCoeffCloseFormQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector);

int main()
{

    // Implementation of a simple minimum snap trajectory generation
    // with polynomial representation
    // using Quadratic Programming with only equality constraints to go through the waypoints

    // 1. takes in waypoints
    // vector<Vector3d> waypt_list{
    //     Vector3d(-2.0, 4.0, 1.0),
    //     Vector3d(-1.0, 0.5, 2.0),
    //     Vector3d(0.0, 0.0, 3.0),
    //     Vector3d(1.0, 1.0, 1.0)};

    log("number of waypts");
    int waypt_num;
    cin >> waypt_num;
    vector<Vector3d> waypt_list(waypt_num);
    
    for(int i = 0; i< waypt_num;i++)
    {
        std::cout << "Enter waypts: "<< (i+1) << std::endl;
        Vector3d pt;
        cin >> pt.x();
        cin >> pt.y();
        cin >> pt.z();
        waypt_list[i] = pt;
    }

    // 3. Set up optimization problem
    // 3.0 get polynomial order (n), derivative order (r) and segment number (k)
    int k_segment = waypt_list.size() - 1;
    int n_poly_order = 5;
    int r_derivate_order = 4; // which derivative's squared integral to minimize

    vector<vector<double>> waypt_single_axis(3);
    for (int j = 0; j < waypt_list.size(); j++)
    {
        waypt_single_axis[0].push_back(waypt_list[j].x());
        waypt_single_axis[1].push_back(waypt_list[j].y());
        waypt_single_axis[2].push_back(waypt_list[j].z());
    }

    // 2. Compute Time allocation for each segment

    double vel_x, vel_y, vel_z;

    vel_x = 2.0;
    vel_y = 2.0;
    vel_z = 2.0;

    double v_norm = 2;

    Vector3d vel_max = {vel_x, vel_y, vel_z};

    vector<double> segment_length(k_segment, 0);

    double distance = computeDistance(waypt_list, segment_length);

    double time = distance/v_norm;

    vector<double> time_vector = computeT(waypt_list, time);

    // vector<double> time_vector = segment_length/distance * time;

    for (auto time : time_vector)
    {
        log(time);
    }

    for (auto length : segment_length)
    {
        log(length);
    }

    // 3.1 decision variables/ vectors that stores the polynomial coefficients
    int coeff_number = k_segment * (n_poly_order + 1);
    VectorXd x_coeff(coeff_number);
    VectorXd y_coeff(coeff_number);
    VectorXd z_coeff(coeff_number);

    // 3.2 objective function

    // f(x) = x' * Q * x

    // SparseMatrix hessian matrix hessian_Q SparseMatrix<double> sp(r_dim, c_dim)
    int hessian_dim = (n_poly_order + 1) * k_segment;
    int decision_dim = hessian_dim;
    // MatrixXd hessian_Q(hessian_dim, hessian_dim, 0);
    // or
    // MatrixXd hessian_Q = MatrixXd::Constant(hessian_dim, hessian_dim, 0.0);
    // or
    // MatrixXd hessian_Q = MatrixXd::Zero(hessian_dim, hessian_dim);
    // or
    // MatrixXd hessian_Q(hessian_dim, hessian_dim);
    MatrixXd hessian_Q = MatrixXd::Zero(hessian_dim, hessian_dim);

    computeHessian(hessian_Q, hessian_dim, n_poly_order, r_derivate_order, k_segment, time_vector);

    // computeCoeffConstrainedQP(x_coeff, hessian_Q, n_poly_order, k_segment, waypt_single_axis[0], time_vector);
    // computeCoeffConstrainedQP(y_coeff, hessian_Q, n_poly_order, k_segment, waypt_single_axis[1], time_vector);
    // computeCoeffConstrainedQP(z_coeff, hessian_Q, n_poly_order, k_segment, waypt_single_axis[2], time_vector);

    computeCoeffCloseFormQP(x_coeff, hessian_Q, n_poly_order, k_segment, waypt_single_axis[0], time_vector);
    computeCoeffCloseFormQP(y_coeff, hessian_Q, n_poly_order, k_segment, waypt_single_axis[1], time_vector);
    computeCoeffCloseFormQP(z_coeff, hessian_Q, n_poly_order, k_segment, waypt_single_axis[2], time_vector);
    log("x_coeff is ");
    log(x_coeff.transpose());
    log("y_coeff is ");
    log(y_coeff.transpose());
    log("z_coeff is ");
    log(z_coeff.transpose());
}

double computeDistance(const vector<Vector3d> &waypt_list, vector<double> &segment_length)
{
    double total_dist = 0;
    for (int i = 0; i < waypt_list.size()-1; i++)
    {
        double length = (waypt_list[i + 1] - waypt_list[i]).norm();
        segment_length[i] = length;
        total_dist += length;
    }

    return total_dist;
}

vector<double> computeT(const vector<Vector3d> &waypt_list, const Vector3d &vel_max)
{
    double vel_norm = vel_max.norm();

    int segment_number = waypt_list.size() - 1;
    int timestamp_number = segment_number + 1;

    vector<double> segment_length(segment_number, 0);
    vector<double> time_vector(segment_number, 0);

    double total_dist = 0;
    for (int i = 0; i < segment_number; i++)
    {
        double length = (waypt_list[i + 1] - waypt_list[i]).norm();
        segment_length[i] = length;
        total_dist += length;
    }

    double total_time = total_dist / vel_norm;

    for (int i = 0; i < segment_number; i++)
    {
        // holds relative timestamp for each segment duration
        // time_vector[i] = total_time * (segment_length[i-1] / total_dist);

        // holds absolute timestamp for each segment duration
        time_vector[i] = total_time * (segment_length[i] / total_dist);
    }

    return time_vector;
}

vector<double> computeT(const vector<Vector3d> &waypt_list, double total_time)
{
    // double vel_norm = vel_max.norm();

    int segment_number = waypt_list.size() - 1;
    int timestamp_number = segment_number + 1;

    vector<double> segment_length(segment_number, 0);
    vector<double> time_vector(segment_number, 0);

    double total_dist = 0;
    for (int i = 0; i < segment_number; i++)
    {
        double length = (waypt_list[i + 1] - waypt_list[i]).norm();
        segment_length[i] = length;
        log(length);
        total_dist += length;
    }

    log(total_dist);

    // double total_time = total_dist / vel_norm;

    for (int i = 0; i < segment_number; i++)
    {
        // holds relative timestamp for each segment duration
        // time_vector[i] = total_time * (segment_length[i-1] / total_dist);

        // holds absolute timestamp for each segment duration
        time_vector[i] = total_time * (segment_length[i] / total_dist);
    }

    return time_vector;
}

void computeHessian(MatrixXd &hessian, int hessian_dim, int n_poly_order, int r_derivate_order, int k_segment, const vector<double> &time_vector)
{
    int sub_mat_dim = n_poly_order + 1;
    // MatrixXd sub_mat = MatrixXd::Zero(sub_mat_dim, sub_mat_dim);
    MatrixXd sub_mat(sub_mat_dim, sub_mat_dim);

    double t_0 = 0;
    double t_f;

    for (int i = 0; i < k_segment; i++) // every segment has one sub_mat
    {
        t_f = time_vector[i];
        t_0 = 0;
        for (int j = r_derivate_order; j < sub_mat_dim; j++)     // row index of sub_mat non-zero entry
        {                                                        // e.g. jerk(t) is 3^rd order derivative
                                                                 // monomial basis becomes [0 0 0 6 24t] for 4^th order poly
                                                                 // non-zero entry starts from index [3]
            for (int k = r_derivate_order; k < sub_mat_dim; k++) // col index of sub_mat non-zero entry
            {
                // before integration
                // sub_mat(j, k) = double(factorial(j)) / factorial(j - r_derivate_order) * pow(t, j - r_derivate_order) *
                //                 double(factorial(k)) / factorial(k - r_derivate_order) * pow(t, k - r_derivate_order);

                // after integration
                sub_mat(j, k) = double(factorial(j)) / factorial(j - r_derivate_order) * double(factorial(k)) / factorial(k - r_derivate_order) / (j - r_derivate_order + k - r_derivate_order + 1) * pow(t_f, j - r_derivate_order + k - r_derivate_order + 1) -
                                double(factorial(j)) / factorial(j - r_derivate_order) * double(factorial(k)) / factorial(k - r_derivate_order) / (j - r_derivate_order + k - r_derivate_order + 1) * pow(t_0, j - r_derivate_order + k - r_derivate_order + 1);
            }
        }

        hessian.block(i * sub_mat_dim, i * sub_mat_dim, sub_mat_dim, sub_mat_dim) = sub_mat;
    }
}

void computeEquality(MatrixXd &Aeq, int Aeq_row, int Aeq_col, VectorXd &beq, int beq_dim, int n_poly_order, int k_segment, const vector<double> &waypoint_list, const vector<double> &time_vector)
{
    // (1) waypoints constraints (k+1)
    // Aeq(0,0) = 1;
    // for (int i = 1; i < k_segment; i++)
    // {
    //     VectorXd vect(n_poly_order + 1);
    //     evaluatePoly(vect, 0, time_vector[i-1], n_poly_order);
    //     Aeq.block(i, i*(n_poly_order + 1), 1, vect.size()) = vect.transpose();
    // }

    // (1) initial p,v,a. row 1 to row 3
    // i = 0, 1, 2
    int row_count = 0;
    beq(row_count) = waypoint_list[0];
    for (int i = 0; i < 3; i++)
    {
        VectorXd vect(n_poly_order + 1);
        evaluatePoly(vect, i, 0.0, n_poly_order);
        Aeq.block(row_count, 0, 1, vect.size()) = vect.transpose();
        row_count++;
    }

    // final p,v,a row 4, 5, 6
    beq(row_count) = waypoint_list.back();
    for (int i = 0; i < 3; i++)
    {
        VectorXd vect(n_poly_order + 1);
        evaluatePoly(vect, i, time_vector[k_segment], n_poly_order);
        Aeq.block(row_count, (k_segment - 1) * (n_poly_order + 1), 1, vect.size()) = vect.transpose();
        row_count++;
    }

    // intermediate waypoints row 7 to row 7 + (k_segment - 1) - 1
    for (int i = 0; i < k_segment - 1; i++)
    {
        VectorXd vect(n_poly_order + 1);
        evaluatePoly(vect, 0, time_vector[i + 1], n_poly_order);
        Aeq.block(row_count, (i + 1) * (n_poly_order + 1), 1, vect.size()) = vect.transpose();
        beq(row_count) = waypoint_list[i + 1];
        row_count++;
    }

    // continuous constraints row 7 + (k_segment - 1) to row 7 + (k_segment - 1) + (k_segment - 1) * 3 - 1
    for (int i = 0; i < k_segment - 1; i++) // for each intermediate waypoint
    {
        for (int j = 0; j < 3; j++) // p,v,a continuous
        {
            VectorXd vect_front(n_poly_order + 1);
            VectorXd vect_back(n_poly_order + 1);
            evaluatePoly(vect_front, j, time_vector[i + 1], n_poly_order);
            vect_back = vect_front * (-1);
            VectorXd vect((n_poly_order + 1) * 2);
            vect << vect_front, vect_back;

            Aeq.block(row_count, i * (n_poly_order + 1), 1, vect.size()) = vect.transpose();
            row_count++;
        }
    }
}

void computeAeq(MatrixXd &Aeq, int n_poly_order, int k_segment, const vector<double> &time_vector)
{
    // (1) waypoints constraints (k+1)
    // Aeq(0,0) = 1;
    // for (int i = 1; i < k_segment; i++)
    // {
    //     VectorXd vect(n_poly_order + 1);
    //     evaluatePoly(vect, 0, time_vector[i-1], n_poly_order);
    //     Aeq.block(i, i*(n_poly_order + 1), 1, vect.size()) = vect.transpose();
    // }

    // (1) initial p,v,a. row 1 to row 3
    // i = 0, 1, 2
    int row_count = 0;
    for (int i = 0; i < 3; i++)
    {
        VectorXd vect(n_poly_order + 1);
        evaluatePoly(vect, i, 0.0, n_poly_order);
        Aeq.block(row_count, 0, 1, vect.size()) = vect.transpose();
        row_count++;
    }

    // final p,v,a row 4, 5, 6
    for (int i = 0; i < 3; i++)
    {
        VectorXd vect(n_poly_order + 1);
        evaluatePoly(vect, i, time_vector[k_segment], n_poly_order);
        Aeq.block(row_count, (k_segment - 1) * (n_poly_order + 1), 1, vect.size()) = vect.transpose();
        row_count++;
    }

    // intermediate waypoints row 7 to row 7 + (k_segment - 1) - 1
    for (int i = 0; i < k_segment - 1; i++)
    {
        VectorXd vect(n_poly_order + 1);
        evaluatePoly(vect, 0, time_vector[i + 1], n_poly_order);
        Aeq.block(row_count, (i + 1) * (n_poly_order + 1), 1, vect.size()) = vect.transpose();
        row_count++;
    }

    // continuous constraints row 7 + (k_segment - 1) to row 7 + (k_segment - 1) + (k_segment - 1) * 3 - 1
    for (int i = 0; i < k_segment - 1; i++) // for each intermediate waypoint
    {
        for (int j = 0; j < 3; j++) // p,v,a continuous
        {
            VectorXd vect_front(n_poly_order + 1);
            VectorXd vect_back(n_poly_order + 1);
            evaluatePoly(vect_front, j, time_vector[i + 1], n_poly_order);
            vect_back = vect_front * (-1);
            VectorXd vect((n_poly_order + 1) * 2);
            vect << vect_front, vect_back;

            Aeq.block(row_count, i * (n_poly_order + 1), 1, vect.size()) = vect.transpose();
            row_count++;
        }
    }
}

void computeBeq(VectorXd &beq, int k_segment, const vector<double> &waypoint_list)
{
    int row_count = 0;
    // [waypoint_list[0], 0, 0, waypoint_list[k_segment], 0, 0, waypoint_list[1], waypoint_list[2], waypoint_list[k_segment-1]]
    beq(row_count) = waypoint_list[0];
    beq(3) = waypoint_list.back();

    for (int i = 0; i < k_segment - 1; i++)
    {
        beq(i + 6) = waypoint_list[i + 1];
    }
}

void evaluatePoly(VectorXd &row_vector, int derivative_order, double time, int polynomial_order)
{
    switch (derivative_order)
    {
    case 0:
    { /* code */
        for (int i = 0; i < polynomial_order + 1; i++)
        {
            row_vector(i) = pow(time, i);
        }
        break;
    }

    case 1:
    { /* code */
        for (int i = 1; i < polynomial_order + 1; i++)
        {
            row_vector(i) = i * pow(time, i - 1);
        }
        break;
    }

    case 2:
    { /* code */
        for (int i = 2; i < polynomial_order + 1; i++)
        {
            row_vector(i) = i * (i - 1) * pow(time, i - 2);
        }
        break;
    }
    default:
        break;
    }
}

VectorXd optimize(const MatrixXd &hessian_Q, const MatrixXd &equality_A, VectorXd &equality_b)
// void optimize(const MatrixXd &hessian_Q, const MatrixXd &equality_A, VectorXd &equality_b)
{
    // int beq_dim = Aeq_row;
    // VectorXd equality_b = VectorXd::Zero(beq_dim);
    // computeBeq(equality_b, beq_dim, k_segment, waypt_list);

    // objective_matrix (SparseMatrix): hessian_Q
    // objective vector (VectorXd): zero()
    // constraint_matrix (SparseMatrix): equality_A
    // lb = ub (VectorXd): equality_b
    // solver.setMats(objective_matrix, objective_vector, constraint_matrix, lower_bound, upper_bound, 1e-3, 1e-3);

    VectorXd objective_vector = VectorXd::Zero(equality_A.cols());

    SparseMatrix<double> objective_matrix = hessian_Q.sparseView();
    SparseMatrix<double> constraint_matrix = equality_A.sparseView();

    IOSQP solver;
    c_int flag = solver.setMats(objective_matrix, objective_vector, constraint_matrix, equality_b, equality_b, 1e-3, 1e-3);
    // c_int flag = solver.setMats(hessian_Q.sparseView(), objective_vector, equality_A.sparseView(), equality_b, equality_b, 1e-3, 1e-3);
    if (flag != 0)
    {
        log("problem non-convex");
    }
    solver.solve();
    return solver.getPrimalSol();
}

void computeDerivateMappingMatrix(MatrixXd &mapping_A, int k_segment, const vector<double> &time_vector)
{
    int sub_dim = mapping_A.cols() / k_segment;

    VectorXd sub_A_row(sub_dim);

    for (int i = 0; i < k_segment; i++)
    {

        MatrixXd sub_A = MatrixXd::Zero(sub_dim, sub_dim);

        for (int j = 0; j < 3; j++)
        {
            evaluatePoly(sub_A_row, j, 0, (sub_dim - 1));
            sub_A.row(j) = sub_A_row.transpose();
            sub_A_row = VectorXd::Zero(sub_dim);
        }

        for (int j = 3; j < 6; j++)
        {
            evaluatePoly(sub_A_row, j - 3, time_vector[i], (sub_dim - 1));
            sub_A.row(j) = sub_A_row.transpose();
            sub_A_row = VectorXd::Zero(sub_dim);
        }
        int block_index = sub_dim * i;

        mapping_A.block(block_index, block_index, sub_dim, sub_dim) = sub_A;
    }
}

void computeCoeffConstrainedQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector)
{
    // // SparseMatrix equality constraint equality_A
    int Aeq_row, Aeq_col;

    // number of rows of Aeq is same as number of equality constraints
    // (1) intermediate waypoints: k-1
    // (2) continuous p,v,a at intermediate waypoints: 3*(k-1)
    // (3) start and end p,v,a: 6
    // total = 4(k-1) + 6 = 4k+2
    // each segment is represented by a n^th order polynomial which has (n+1) coefficients
    Aeq_row = 4 * k_segment + 2;

    // number of columns of Aeq is same as number of decision variables (number of polynomial coefficients)
    // each segment is represented by a n^th order polynomial which has (n+1) coefficients
    Aeq_col = k_segment * (n_poly_order + 1);

    MatrixXd equality_A = MatrixXd::Zero(Aeq_row, Aeq_col);
    int beq_dim = Aeq_row;
    VectorXd equality_b = VectorXd::Zero(beq_dim);

    computeAeq(equality_A, n_poly_order, k_segment, time_vector);

    computeBeq(equality_b, k_segment, waypt_single_axis);

    coeff = optimize(hessian_Q, equality_A, equality_b);
}

void computeCoeffCloseFormQP(VectorXd &coeff, const MatrixXd &hessian_Q, int n_poly_order, int k_segment, const vector<double> &waypt_single_axis, const vector<double> &time_vector)
{
    // closed form solutions
    int mapping_A_dim = k_segment * 6; // every segment will have 2 sets of p,v,a constraints
    MatrixXd mapping_A = MatrixXd::Zero(mapping_A_dim, mapping_A_dim);
    computeDerivateMappingMatrix(mapping_A, k_segment, time_vector);

    //
    int constraint_order = 2;                                      // p (order 0), v (order 1), a (order 2)
    int fixed_number = 2 * (constraint_order + 1) + k_segment - 1; // initial and final p,v,a + intermediate p
    int free_number = constraint_order * (k_segment - 1);          // intermediate v,a

    // VectorXd decision_f(fixed_number), decision_p(free_number);

    VectorXd decision_f = MatrixXd::Zero(fixed_number, 1);
    decision_f(0) = waypt_single_axis[0];
    for (int i = 3; i < 3 + k_segment; i++)
    {
        decision_f(i) = waypt_single_axis[i - 2];
    }

    VectorXd decision_p_primal = MatrixXd::Zero(free_number, 1);

    int C_t_col = fixed_number + free_number;
    int C_t_row = k_segment * 2 * (1 + constraint_order);

    MatrixXd C_t = MatrixXd::Zero(C_t_row, C_t_col);

    Matrix3d C_0;
    C_0 << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;

    Vector4d C_1;
    C_1 << 1, 0, 0, 1;

    Matrix<double, 5, 2> C_2;
    C_2.setZero();
    C_2(0, 0) = 1;
    C_2(1, 1) = 1;
    C_2(3, 0) = 1;
    C_2(4, 1) = 1;

    int row_count = 0;
    int col_count = 0;

    C_t.block(row_count, col_count, C_0.rows(), C_0.cols()) = C_0;

    row_count += C_0.rows();
    col_count += C_0.cols();

    for (int i = 0; i < k_segment - 1; i++)
    {
        C_t.block(row_count, col_count, C_1.rows(), C_1.cols()) = C_1;
        row_count += C_1.rows();
        row_count += 2;
        col_count += C_1.cols();
    }

    C_t.block(row_count, col_count, C_0.rows(), C_0.cols()) = C_0;
    col_count += C_0.cols();
    row_count = 4;

    for (int i = 0; i < k_segment - 1; i++)
    {
        C_t.block(row_count, col_count, C_2.rows(), C_2.cols()) = C_2;
        row_count += C_2.rows();
        row_count += 1;
        col_count += C_2.cols();
    }

    int dim_R = fixed_number + free_number;
    MatrixXd R = MatrixXd::Zero(dim_R, dim_R);
    R = C_t.transpose() * (mapping_A.inverse()).transpose() * hessian_Q * mapping_A.inverse() * C_t;

    MatrixXd R_ff = MatrixXd::Zero(fixed_number, fixed_number);
    MatrixXd R_fp = MatrixXd::Zero(fixed_number, free_number);
    MatrixXd R_pf = MatrixXd::Zero(free_number, fixed_number);
    MatrixXd R_pp = MatrixXd::Zero(free_number, free_number);

    R_ff = R.block(0, 0, fixed_number, fixed_number);
    R_fp = R.block(0, fixed_number, fixed_number, free_number);
    R_pf = R.block(fixed_number, 0, free_number, fixed_number);
    R_pp = R.block(fixed_number, fixed_number, free_number, free_number);

    decision_p_primal = -R_pp.inverse() * R_fp.transpose() * decision_f;
    VectorXd decision_primal(fixed_number + free_number);
    decision_primal << decision_f, decision_p_primal;

    VectorXd derivative_vector(k_segment * (constraint_order + 1) * 2);

    derivative_vector = C_t * decision_primal;

    coeff = mapping_A.inverse() * derivative_vector;
}