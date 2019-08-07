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
    bool SolvePnPByEPNP(const std::vector<Eigen::Vector3d> &pts3d,const std::vector<Eigen::Vector2d> &pts2d,const Eigen::Matrix3d & K,Eigen::Matrix3d &R,Eigen::Vector3d & t);



private:

    //---These functions and variables are used by DLT
    bool ConstructM(const std::vector<Eigen::Vector3d >& pts3d, const std::vector< Eigen::Vector2d >& pts2d);
    bool ComputeP();
    bool DecomposeP(Eigen::Matrix3d& R, Eigen::Vector3d& t, Eigen::Matrix3d& K);
    Eigen::Matrix<double,3,4> P_;
    Eigen::MatrixXd M_;
    //-----------------



    //---The following  Functions or variables are used by EPNP
    void ChooseControlPointsWorld(const std::vector<Eigen::Vector3d >& pts3d,std::vector<Eigen::Vector3d> &control_points_world);
    void ComputeHBCoordinates(const std::vector<Eigen::Vector3d> & pts3d, const std::vector<Eigen::Vector3d> &control_points_world , std::vector<Eigen::Vector4d> &hb_coordinates);
    void ComputeM(const std::vector<Eigen::Vector2d>& pts2d, const Eigen::Matrix3d&  K, const std::vector<Eigen::Vector4d>& hbcs, Eigen::MatrixXd &M  );
    void ComputeL6x10(const Eigen::Matrix<double ,12,4> &  V,Eigen::Matrix<double,6,10> & L6x10);
    void ComputeRho(const std::vector<Eigen::Vector3d> &control_points_world ,Eigen::VectorXd &rho );
    void ComputeBetaApproximate1(const Eigen::Matrix<double ,6,10> & L6x10,const Eigen::VectorXd& rho,Eigen::Vector4d& beta);
    void DoGaussNewtonOptimization(const Eigen::Matrix<double ,6,10>&  L6x10,const Eigen::VectorXd& rho,Eigen::Vector4d& beta);
    void ComputeJacobianMatrix(const Eigen::Matrix<double ,6,10>&  L6x10,const Eigen::Vector4d& beta ,Eigen::Matrix<double ,6,4> & jacobian);
    void ComputeResidual(const Eigen::Matrix<double ,6,10>&  L6x10,const Eigen::Vector4d& beta,const Eigen::VectorXd& rho,Eigen::VectorXd& residual );

    void ComputeControlPointsCamera(const Eigen::Vector4d& beta,const Eigen::Matrix<double ,12,4> &V,std::vector<Eigen::Vector3d >& control_points_camera);
    void ComputePts3DCameraCoordinates(const std::vector<Eigen::Vector4d>& hb_coordinates, const std::vector<Eigen::Vector3d >& control_points_camera, std::vector<Eigen::Vector3d>& pts_3d_camera );
    void ComputeRt(const std::vector<Eigen::Vector3d>& pts3d,const std::vector<Eigen::Vector3d > & pts_3d_camera,Eigen::Matrix3d& R,Eigen::Vector3d& t);





};



#endif //PNP_SOLVER_PNP_SOLVER_H
