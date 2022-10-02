#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#define log(x) std::cout << x << std::endl;

using namespace std;
using namespace Eigen;

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

int main()
{

    int derivative_order = 0;
    double time = 0;
    int polynomial_order = 5;
    MatrixXd mat(10, 10);

    log(mat);
    VectorXd vect(polynomial_order + 1);
    evaluatePoly(vect, derivative_order, time, polynomial_order);
    log("vect is ")
        log(vect);

    mat.block(0, 0, 1, vect.size()) = vect.transpose();
    vect = VectorXd::Zero(vect.size());
    evaluatePoly(vect, 1, time, polynomial_order);
    log("vect is ")
        log(vect);

    mat.block(1, 0, 1, vect.size()) = vect.transpose();
    vect = VectorXd::Zero(vect.size());
    evaluatePoly(vect, 2, time, polynomial_order);
    log("vect is ")
        log(vect);
    mat.block(2, 0, 1, vect.size()) = vect.transpose();
    log(mat);

    VectorXd vect1(5);
    vect1 << 1, 2, 3, 4, 5;

    log(vect1 * 2);

    VectorXd vect2(10);
    vect2 << vect1, -vect1;
    log(vect2.transpose());

    vector<Vector3d> waypt_list{
        Vector3d(-2.0, 4.0, 1.0),
        Vector3d(-1.0, 0.5, 1.0),
        Vector3d(0.0, 0.0, 1.0),
        Vector3d(1.0, 1.0, 1.0)};

    vector<vector<double>> waypt_single_axis(3);

    for (int j = 0; j < waypt_list.size(); j++)
    {
        waypt_single_axis[0].push_back(waypt_list[j].x());
        waypt_single_axis[1].push_back(waypt_list[j].y());
        waypt_single_axis[2].push_back(waypt_list[j].z());
    }

    for (auto & pt : waypt_single_axis[0])
    {
        log(pt);
    }

        for (auto & pt : waypt_single_axis[1])
    {
        log(pt);
    }

        for (auto & pt : waypt_single_axis[2])
    {
        log(pt);
    }
}