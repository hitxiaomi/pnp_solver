//
// Created by xiaomi on 8/1/19.
//

#ifndef PNP_SOLVER_PNP_SOLVER_H
#define PNP_SOLVER_PNP_SOLVER_H
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
class PnPSolver{
public:
    bool SolvePnPByDLT( const std::vector< Eigen::Vector3d >& pts3d, const std::vector< Eigen::Vector2d >& pts2d, Eigen::Matrix3d& R, Eigen::Vector3d& t, Eigen::Matrix3d& K);




private:
    bool ConstructM(const std::vector<Eigen::Vector3d >& pts3d, const std::vector< Eigen::Vector2d >& pts2d);
    bool ComputeP();
    bool DecomposeP(Eigen::Matrix3d& R, Eigen::Vector3d& t, Eigen::Matrix3d& K);
    Eigen::Matrix<double,3,4> P_;
    Eigen::MatrixXd M_;
};



#endif //PNP_SOLVER_PNP_SOLVER_H
