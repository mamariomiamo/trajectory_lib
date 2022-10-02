#include "matplotlibcpp.h"
#include "iosqp.hpp"
#include <iostream>
#include <eigen3/Eigen/Dense>
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

int main()
{

        vector<Vector3d> waypt_list{
            Vector3d(-2.0, 4.0, 1.0),
            Vector3d(-1.0, 0.5, 1.0),
            Vector3d(0.0, 0.0, 1.0),
            Vector3d(1.0, 1.0, 1.0)};

        waypt_list.push_back(Vector3d(0.0, 0.0, 0.0));

        std::cout << waypt_list[0] << std::endl;
        std::cout << waypt_list.size() << std::endl;

    //     // plt::plot({1,3,2,4});
    //     // plt::show();

    //     // plt::scatter(waypt_list[0].x(), waypt_list[0].y(); waypt_list[0].z());
    //     // plt::show();

    //     std::vector<std::vector<double>> x, y, z;
    //     for (double i = -5; i <= 5; i += 0.25)
    //     {
    //         std::vector<double> x_row, y_row, z_row;
    //         for (double j = -5; j <= 5; j += 0.25)
    //         {
    //             x_row.push_back(i);
    //             y_row.push_back(j);
    //             z_row.push_back(::std::sin(::std::hypot(i, j)));
    //         }
    //         x.push_back(x_row);
    //         y.push_back(y_row);
    //         z.push_back(z_row);
    //     }

    //     // plt::plot_surface(x, y, z);
    //     // plt::show();

        vector<double> waypt_x, waypt_y, waypt_z;
        for (int i = 0; i < waypt_list.size(); i++)
        {
            waypt_x.push_back(waypt_list[i].x());
            waypt_y.push_back(waypt_list[i].y());
            waypt_z.push_back(waypt_list[i].z());
        }

        plt::plot3(waypt_x, waypt_y, waypt_z);
        plt::show();

        plt::plot(waypt_x, waypt_y);
        plt::show();

        plt::plot(waypt_x);
        plt::show();

        double nominal_vel = 5;
        Vector3d max_vel{3.0, 3.0, 3.0};

        vector<vector<double>> poly_coeff(3, {1.3, 0.0, 0.0, 0.0, 0.0});
        // poly_coeff.reserve(3);

        // for(int i=0; i <3; i++)
        // {
        //     log(poly_coeff[i][0]);
        // }

        // log(poly_coeff.capacity());

        Matrix<double, 6, 6> q_hat = MatrixXd::Constant(6, 6, 0.0);

        Matrix<double, 12, 12> q = MatrixXd::Constant(12, 12, 0.0);
        // Matrix<double, 5, 5> q_hat;
        // q_hat.row(0) << 1,2,3,4,5;
        // q_hat(0) = {1, 2, 3, 4, 5};
        // q_hat = MatrixXd::Constant(5,5,0.0);
        int r_derivative = 4;
        int n_order = 5;
        int n_segments;

        for (int i = r_derivative; i < n_order + 1; i++)
        {
            for (int j = r_derivative; j < n_order + 1; j++)
            {
                q_hat(i, j) = double(factorial(i)) / factorial(i - r_derivative) * pow(1, i - r_derivative) *
                              double(factorial(j)) / factorial(j - r_derivative) * pow(1, j - r_derivative);
            }
        }

        q.block<6, 6>(0, 0) = q_hat;
        q.block<6, 6>(6, 6) = q_hat;

        log(q);

    int m_derivative;
    int m_order;
    int m_segments;

    log("derivative order?");
    std::cin >> m_derivative;

    log("polynomial order?");
    std::cin >> m_order;

    log("number of segments?");
    std::cin >> m_segments;

    int m_dim = m_order + 1;
    int m_q_dim = m_dim * m_segments;

    MatrixXd small_q(m_dim, m_dim);
    MatrixXd hessian_q(m_q_dim, m_q_dim);

    for (int k = 0; k < m_segments; k++)
    {
        for (int i = m_derivative; i < m_dim; i++)
        {
            for (int j = m_derivative; j < m_dim; j++)
            {
                small_q(i, j) = double(factorial(i)) / factorial(i - m_derivative) * pow(1, i - m_derivative) *
                                double(factorial(j)) / factorial(j - m_derivative) * pow(1, j - m_derivative);
            }
        }

        // hessian_q.topLeftCorner<m_dim, m_dim>(k * m_dim) = small_q;
        hessian_q.block(k * m_dim, k * m_dim, m_dim, m_dim) = small_q;
        log(small_q);
    }

    // log(small_q);
    // log(hessian_q);

    IOSQP solver;
    SparseMatrix<double> P(2, 2);
    const Triplet<double> kTripletsP[] = {
        {0, 0, 4.0}, {1, 0, 1.0}, {0, 1, 1.0}, {1, 1, 2.0}};
    P.setFromTriplets(std::begin(kTripletsP),
                      std::end(kTripletsP));

    log(P);

    // VectorXd q_(1.0, 1.0);
    Vector3d a;
    a << 1.0, 1.0;

    SparseMatrix<double> A(3, 2);

    // log(A);

    // log(q_);

    // testQuadProg();
}

// int main()
// {
//     int r_derivative;
//     int n_order;
//     int n_segments;

//     log("derivative order?");
//     std::cin >> r_derivative;

//     log("polynomial order?");
//     std::cin >> n_order;

//     log("number of segments?");
//     std::cin >> n_segments;

//     int matrix_size = n_order + 1;
//     int q_dim = matrix_size * n_segments;

//     MatrixXd<double, matrix_size, matrix_size> q_hat = MatrixXd::Constant(matrix_size, matrix_size, 0.0);

//     MatrixXd<double, q_dim, q_dim> q = MatrixXd::Constant(q_dim, q_dim, 0.0);
//     // Matrix<double, 5, 5> q_hat;
//     // q_hat.row(0) << 1,2,3,4,5;
//     // q_hat(0) = {1, 2, 3, 4, 5};
//     // q_hat = MatrixXd::Constant(5,5,0.0);

//     for (k = 0; k < n_segments; k++)
//     {
//         for (int i = r_derivative; i < matrix_size; i++)
//         {
//             for (int j = r_derivative; j < matrix_size; j++)
//             {
//                 q_hat(i, j) = double(factorial(i)) / factorial(i - r_derivative) * pow(1, i - r_derivative) *
//                               double(factorial(j)) / factorial(j - r_derivative) * pow(1, j - r_derivative);
//             }
//         }

//         q.block<matrix_size, matrix_size>(k * matrix_size) = q_hat;
//     }

//     // q.block<6, 6>(0, 0) = q_hat;
//     // q.block<6, 6>(6, 6) = q_hat;

//     log(q);
// }