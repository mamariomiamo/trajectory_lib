#include "iosqp.hpp"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <iterator>

using namespace std;
using namespace Eigen;

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
    Eigen::SparseMatrix<double> objective_matrix(2, 2);
    const Triplet<double> kTripletsP[] = {
        {0, 0, 1}, {1, 0, -1}, {0, 1, -1}, {1, 1, 2}};
    objective_matrix.setFromTriplets(std::begin(kTripletsP),
                                     std::end(kTripletsP));

    Eigen::SparseMatrix<double> constraint_matrix(3, 2);
    const Triplet<double> kTripletsA[] = {
        {0, 0, 1}, {1, 0, -1}, {2, 0, 2}, {0, 1, 1}, {1, 1, 2}, {2, 1, 1}};
    constraint_matrix.setFromTriplets(std::begin(kTripletsA),
                                      std::end(kTripletsA));

    std::cout << "constraint is\n"
              << constraint_matrix << std::endl;
    Eigen::VectorXd objective_vector(2);
    objective_vector << -2, -6;

    Eigen::VectorXd upper_bound(3);
    upper_bound << 2.0, 2.0, 3.0;

    Eigen::VectorXd lower_bound(3);
    lower_bound << -kInfinity, -kInfinity, -kInfinity;

    IOSQP solver(false);
    c_int flag = solver.setMats(objective_matrix, objective_vector, constraint_matrix, lower_bound, upper_bound, 1e-3, 1e-3);
    if (flag != 0)
    {
        std::cout << "problem non-convex\n";
    }
    solver.solve();
    // c_int status = solver.getStatus();
    //  std::cout << "STATUS: " << status << std::endl;
    // q_ << 1.0, 1.0;

    std::cout << "primal solution is " << solver.getPrimalSol() << std::endl;
}

int main()
{

    std::vector<Eigen::Vector3d> waypt_list{
        Eigen::Vector3d(-2.0, 4.0, 1.0),
        Eigen::Vector3d(-1.0, 0.5, 1.0),
        Eigen::Vector3d(0.0, 0.0, 1.0),
        Eigen::Vector3d(1.0, 1.0, 1.0)};

    waypt_list.push_back(Eigen::Vector3d(0.0, 0.0, 0.0));

    std::cout << waypt_list[0] << std::endl;
    std::cout << waypt_list.size() << std::endl;

    std::vector<double> waypt_x, waypt_y, waypt_z;
    for (int i = 0; i < waypt_list.size(); i++)
    {
        waypt_x.push_back(waypt_list[i].x());
        waypt_y.push_back(waypt_list[i].y());
        waypt_z.push_back(waypt_list[i].z());
    }

    double nominal_vel = 5;
    Eigen::Vector3d max_vel{3.0, 3.0, 3.0};

    std::vector<vector<double>> poly_coeff(3, {1.3, 0.0, 0.0, 0.0, 0.0});

    Eigen::Matrix<double, 6, 6> q_hat = Eigen::MatrixXd::Constant(6, 6, 0.0);

    Eigen::Matrix<double, 12, 12> q = Eigen::MatrixXd::Constant(12, 12, 0.0);

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

    std::cout << "q is " << q << std::endl;

    int m_derivative;
    int m_order;
    int m_segments;

    std::cout << "derivative order?\n";
    std::cin >> m_derivative;

    std::cout <<"polynomial order?\n";
    std::cin >> m_order;

    std::cout <<"number of segments?\n";
    std::cin >> m_segments;

    int m_dim = m_order + 1;
    int m_q_dim = m_dim * m_segments;

    Eigen::MatrixXd small_q(m_dim, m_dim);
    Eigen::MatrixXd hessian_q(m_q_dim, m_q_dim);

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
        std::cout << small_q << std::endl;
    }

    IOSQP solver(false);
    Eigen::SparseMatrix<double> P(2, 2);
    const Triplet<double> kTripletsP[] = {
        {0, 0, 4.0}, {1, 0, 1.0}, {0, 1, 1.0}, {1, 1, 2.0}};
    P.setFromTriplets(std::begin(kTripletsP),
                      std::end(kTripletsP));

    std::cout << P << std::endl;

    // Eigen::VectorXd q_(1.0, 1.0);
    Eigen::Vector3d a;
    a << 1.0, 1.0;

    Eigen::SparseMatrix<double> A(3, 2);

    // std::cout <<A);

    // std::cout <<q_);

    // testQuadProg();
}

// int main()
// {
//     int r_derivative;
//     int n_order;
//     int n_segments;

//     std::cout <<"derivative order?");
//     std::cin >> r_derivative;

//     std::cout <<"polynomial order?");
//     std::cin >> n_order;

//     std::cout <<"number of segments?");
//     std::cin >> n_segments;

//     int matrix_size = n_order + 1;
//     int q_dim = matrix_size * n_segments;

//     Eigen::MatrixXd<double, matrix_size, matrix_size> q_hat = Eigen::MatrixXd::Constant(matrix_size, matrix_size, 0.0);

//     Eigen::MatrixXd<double, q_dim, q_dim> q = Eigen::MatrixXd::Constant(q_dim, q_dim, 0.0);
//     // Eigen::Matrix<double, 5, 5> q_hat;
//     // q_hat.row(0) << 1,2,3,4,5;
//     // q_hat(0) = {1, 2, 3, 4, 5};
//     // q_hat = Eigen::MatrixXd::Constant(5,5,0.0);

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

//     std::cout <<q);
// }