#include "matplotlibcpp.h"
#include "min_snap.h"
#include "iosqp.hpp"
#include <iostream>
// #include <eigen3/Eigen/Dense>
// #include <eigen3/Eigen/Sparse>
#include <vector>
#include <iterator>
#include <memory>
#include <chrono>
#define log(x) std::cout << x << std::endl;

using namespace std;
using namespace Eigen;
using namespace std::chrono;
namespace plt = matplotlibcpp;
// typedef std::chrono::time_point<std::chrono::system_clock> sys_clock;
int main()
{
    // sys_clock now = system_clock::now();
    // std::cout
    //     << "Slow calculations took "
    //     << std::chrono::duration_cast<std::chrono::microseconds> now.count();
    // log(duration_cast<microseconds>now.count());

    vector<Vector3d> waypt_list{
        Vector3d(-2.0, 4.0, 1.0),
        Vector3d(-1.0, 0.5, 2.0),
        Vector3d(0.0, 0.0, 3.0),
        Vector3d(1.0, 1.0, 1.0)};

    const auto close_form_tic = high_resolution_clock::now();
    auto min_snap_close_form = std::make_shared<MIN_SNAP::min_snap_traj>(waypt_list, 3, true);
    // std::shared_ptr<MIN_SNAP::min_snap_traj> min_snap_close_form = std::make_shared<MIN_SNAP::min_snap_traj>(waypt_list, 3, true);
    MIN_SNAP::poly_traj trajectory;
    trajectory = min_snap_close_form->computeTrajectory(waypt_list);
    const auto close_form_toc = high_resolution_clock::now();
    auto close_form_time_lapsed = duration<double>(close_form_toc - close_form_tic).count();
    log("close form computation time in (ms):");
    log(close_form_time_lapsed*1000);
    for (auto entry : trajectory.coeff)
    {
        log(entry);
        log("\n");
    }

    const auto qp_tic = high_resolution_clock::now();
    auto min_snap_qp = std::make_shared<MIN_SNAP::min_snap_traj>(waypt_list, 3, false);
    MIN_SNAP::poly_traj trajectory_qp;
    trajectory_qp = min_snap_qp->computeTrajectory(waypt_list);
    const auto qp_toc = high_resolution_clock::now();
    auto qp_time_lapsed = duration<double>(qp_toc - qp_tic).count();

    log("QP computation time in (ms):");
    log(qp_time_lapsed*1000);
    for (auto entry : trajectory_qp.coeff)
    {
        log(entry);
        log("\n");
    }
}