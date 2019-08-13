#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "pnp_solver.h"
#include <fstream>
#include <random>
bool GenerateTestData(const Eigen::Matrix3d & R,const Eigen::Vector3d & t, const Eigen::Matrix3d &K,std::vector<Eigen::Vector3d> &pts_3d, std::vector<Eigen::Vector2d>  &pts_2d ){

    double w_sigma= 1.0;                 // 噪声Sigma值
    std::default_random_engine generator;
    std::normal_distribution<double> noise(0.,w_sigma);


    std::ifstream infile("3d_points.txt");

    double x,y,z;
    while(infile>>x>>y>>z){
    pts_3d.emplace_back(x,y,z);
    }

    for(const auto &pt:pts_3d){

        double n=noise(generator);
        Eigen::Vector3d tem_p= K*R *(pt-t);
        double u=tem_p(0)/tem_p(2)+n;
        double v=tem_p(1)/tem_p(2)+n;
        pts_2d.emplace_back(u,v);


    }
    return true;


}



int main(int argc,char** argv) {

    Eigen::AngleAxisd rotation_angle_axis(M_PI/4,Eigen::Vector3d(0,0,1));
    Eigen::Matrix3d R_cw=rotation_angle_axis.toRotationMatrix();

    Eigen::Vector3d t_w(1,0,0);
    Eigen::Matrix3d K;
    K<<100,0,160,
        0,100,120 ,
        0,0,1;

    std::vector<Eigen::Vector3d> pts_3d;
    std::vector<Eigen::Vector2d> pts_2d;
    GenerateTestData(R_cw,t_w,K,pts_3d,pts_2d);


    PnPSolver pnp_solver;

    Eigen::Matrix3d R_estimate,K_estimate;
    Eigen::Vector3d t_estimate;

     pnp_solver.SolvePnPByDLT(pts_3d,pts_2d,R_estimate,t_estimate,K_estimate);



     std::cout<<"R"<<std::endl;
     std::cout<<R_cw<<std::endl;

     std::cout<<"R_estimated "<< std::endl;
     std::cout<<R_estimate<<std::endl;
     std::cout<<"R_derter:"<<R_estimate.determinant()<<std::endl;

     std::cout<<"K"<<std::endl;
     std::cout<<K<<std::endl;

     std::cout<<"K_estimated"<<std::endl;
     std::cout<<K_estimate<<std::endl;


    std::cout<<"t_w"<<std::endl;
    std::cout<<t_w<<std::endl;

    std::cout<<"t_estimated"<<std::endl;
    std::cout<<t_estimate<<std::endl;


    return 0;
}