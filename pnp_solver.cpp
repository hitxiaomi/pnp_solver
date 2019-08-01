//
// Created by xiaomi on 8/1/19.
//

#include "pnp_solver.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Geometry>
#include <iostream>
bool PnPSolver::SolvePnPByDLT(const std::vector<Eigen::Vector3d> &pts3d, const std::vector<Eigen::Vector2d> &pts2d,
                              Eigen::Matrix3d &R, Eigen::Vector3d &t, Eigen::Matrix3d &K) {
    if(!ConstructM(pts3d,pts2d))
        return false;
    ComputeP();
    DecomposeP(R,t,K);
    return  true;

}
bool PnPSolver::ConstructM(const std::vector<Eigen::Vector3d> &pts3d, const std::vector<Eigen::Vector2d> &pts2d) {
    int num_of_3d_points=pts3d.size();
    int num_of_2d_points=pts2d.size();
    if((num_of_3d_points < 6) || (num_of_3d_points!= num_of_2d_points)){
        return false;
    }

    M_.resize(2*num_of_3d_points,12);
    for (int i=0;i<num_of_3d_points;i++){
        Eigen::Matrix<double,2,12> matrix_12;
        double x=pts3d[i](0) , y=pts3d[i](1) , z=pts3d[i](2);
        double u=pts2d[i](0) , v=pts2d[i](1);
       matrix_12<<-x,-y,-z,-1,0,0,0,0,u*x,u*y,u*z,u,
                0,0,0,0,-x,-y,-z,-1,v*x,v*y,v*z,v;

        M_.block(2*i,0,2,12)=matrix_12;

    }
    std::cout<<"M"<<std::endl;
    std::cout<<M_<<std::endl;
    return  true;
}
bool PnPSolver::ComputeP() {

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_m(M_,Eigen::ComputeFullV);
    Eigen::MatrixXd V=svd_m.matrixV();


    std::cout<<"---V---:"<<std::endl;
    std::cout<<V<<std::endl;

    Eigen::MatrixXd  Sigma=svd_m.singularValues();
    std::cout<<"---sigma---"<<std::endl;
    std::cout<<Sigma<<std::endl;

    Eigen::Matrix<double ,12,1> p;
    p<<V(0,11),V(1,11),V(2,11),V(3,11),
            V(4,11),V(5,11),V(6,11),V(7,11),
            V(8,11),V(9,11),V(10,11),V(11,11);
    p.normalize();

    P_<<p(0),p(1),p(2),p(3),
    p(4), p(5),p(6),p(7),
    p(8),p(9),p(10),p(11);
    return  true;

}
bool PnPSolver::DecomposeP(Eigen::Matrix3d &R, Eigen::Vector3d &t, Eigen::Matrix3d &K) {

    Eigen::Matrix3d H=P_.block(0,0,3,3);
    Eigen::Matrix3d H_inv=H.inverse();
    Eigen::Vector3d h=P_.col(3);
    Eigen::Vector3d tem_t=-H_inv*h;
    t=tem_t;



    Eigen::HouseholderQR<Eigen::Matrix3d> qr(H_inv);
    Eigen::Matrix3d temR=qr.householderQ();
    Eigen::Matrix3d temK=qr.matrixQR().triangularView<Eigen::Upper>();

    std::cout<<"H_inv:"<<H_inv<<std::endl;
    std::cout<<"Q*R"<<temR*temK<<std::endl;




    R=temR.transpose();
    K=temK.inverse();
    K=K*(1/K(2,2));

    if(K(0,0)<0){

    Eigen::AngleAxisd rotation(M_PI,Eigen::Vector3d::UnitZ());
    Eigen::Matrix3d R_z_pi=rotation.toRotationMatrix() ;
    std::cout<<"R_z_pi"<<std::endl;
    std::cout<<R_z_pi<<std::endl;

    R=R_z_pi*R;
    K=K*R_z_pi;

    }

    return  true;

}