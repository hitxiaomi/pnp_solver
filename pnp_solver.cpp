//
// Created by xiaomi on 8/1/19.
//

#include "pnp_solver.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Geometry>
#include <iostream>
#include <Eigen/Eigenvalues>
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

    return  true;
}
bool PnPSolver::ComputeP() {

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_m(M_,Eigen::ComputeFullV);
    Eigen::MatrixXd V=svd_m.matrixV();
    Eigen::MatrixXd  Sigma=svd_m.singularValues();
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

//    std::cout<<"H_inv:"<<H_inv<<std::endl;
//    std::cout<<"Q*R"<<temR*temK<<std::endl;


    R=temR.transpose();
    K=temK.inverse();
    K=K*(1/K(2,2));

    if(K(0,0)<0){

    Eigen::AngleAxisd rotation(M_PI,Eigen::Vector3d::UnitZ());
    Eigen::Matrix3d R_z_pi=rotation.toRotationMatrix() ;
//    std::cout<<"R_z_pi"<<std::endl;
//    std::cout<<R_z_pi<<std::endl;

    R=R_z_pi*R;
    K=K*R_z_pi;

    }

    return  true;

}
bool PnPSolver::SolvePnPByEPNP(const std::vector<Eigen::Vector3d> &pts3d, const std::vector<Eigen::Vector2d> &pts2d,
                               const Eigen::Matrix3d &K, Eigen::Matrix3d &R, Eigen::Vector3d &t) {

    std::vector<Eigen::Vector3d> control_points_world;

    ChooseControlPointsWorld(pts3d, control_points_world);
//    std::cout<<"control points:"<<std::endl;
//    std::cout<<"----"<<std::endl;
//    for(const auto & pt:control_points_world ){
//        std::cout<<pt<<std::endl;
//    }
//    std::cout<<"----"<<std::endl;



    std::vector<Eigen::Vector4d> hb_coordinates;

    ComputeHBCoordinates(pts3d, control_points_world, hb_coordinates);

    //test hb coordinates and control points
//    std::cout<<"pw combine hb control points---------"<<std::endl;
//    for(const auto & hb: hb_coordinates){
//        Eigen::Vector3d pt_w=hb(0)*control_points_world[0] +hb(1)*control_points_world[1]+hb(2)*control_points_world[2] +hb(3)*control_points_world[3];
//        std::cout<<pt_w<<std::endl<<std::endl;
//    }
    Eigen::MatrixXd M;
    ComputeM(pts2d, K, hb_coordinates, M);
//    std::cout << "M:" << std::endl;
//    std::cout<<M<<std::endl;
    Eigen::Matrix<double, 12, 12> MtM = M.transpose() * M;
//    std::cout << MtM << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 12, 12> > saes(MtM);
    Eigen::MatrixXd eigen_vectors = saes.eigenvectors();
    Eigen::VectorXd eigen_values = saes.eigenvalues();

    std::cout << "eigen_vector:" << std::endl;
    std::cout << eigen_vectors << std::endl;

    Eigen::Matrix<double, 12, 4> V = eigen_vectors.block(0, 0, 12, 4);

    Eigen::Matrix<double, 6, 10> L6x10;
    ComputeL6x10(V, L6x10);
    Eigen::VectorXd rho;
    ComputeRho(control_points_world, rho);

    Eigen::Vector4d beta;

    ComputeBetaApproximate1(L6x10, rho, beta);
    DoGaussNewtonOptimization(L6x10, rho, beta);

    std::vector<Eigen::Vector3d> control_points_camera;
    ComputeControlPointsCamera(beta, V, control_points_camera);

    std::vector<Eigen::Vector3d> pts_3d_camera;
    ComputePts3DCameraCoordinates(hb_coordinates, control_points_camera, pts_3d_camera);

    ComputeRt(pts3d, pts_3d_camera, R, t);

}
void PnPSolver::ChooseControlPointsWorld(const std::vector<Eigen::Vector3d> &pts3d, std::vector<Eigen::Vector3d> &control_points_world) {

    control_points_world.resize(4);
    const int num_of_pts=pts3d.size();
    Eigen::Vector3d sum_of_pts(0,0,0);
    for (const auto pt:pts3d)
        sum_of_pts+=pt;
    control_points_world[0]= sum_of_pts / num_of_pts;

    Eigen::MatrixXd A;
    A.resize(num_of_pts,3);

    for (int i=0;i<num_of_pts;i++){
        A.row(i)= (pts3d[i] - control_points_world[0]).transpose();
    }
    Eigen::Matrix3d AtA= A.transpose() * A;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(AtA);

    Eigen::Matrix3d eigen_vectors=saes.eigenvectors();
    Eigen::Vector3d eigen_values=saes.eigenvalues();

    control_points_world[1]=control_points_world[0] + sqrt(eigen_values(0) /num_of_pts)* eigen_vectors.col(0);
    control_points_world[2]=control_points_world[0] + sqrt(eigen_values(1) /num_of_pts)* eigen_vectors.col(1);
    control_points_world[3]=control_points_world[0] + sqrt(eigen_values(2) /num_of_pts)* eigen_vectors.col(2);



}
void PnPSolver::ComputeHBCoordinates(const std::vector<Eigen::Vector3d> &pts3d,
                                     const std::vector<Eigen::Vector3d> &control_points_world, std::vector<Eigen::Vector4d> &hb_coordinates) {
    Eigen::Matrix4d C=Eigen::Matrix4d::Ones();
    const int num_of_pts3d=pts3d.size();
    hb_coordinates.resize(num_of_pts3d);


    for(int i=0;i<4;i++){
        C.block(0,i,3,1)=control_points_world[i];
    }
    Eigen::Matrix4d C_inv=C.inverse();

    for(int i=0;i<num_of_pts3d;i++){
        Eigen::Vector4d pt_homogeneous=Eigen::Vector4d::Ones();
        pt_homogeneous.head(3)=pts3d[i];
        hb_coordinates[i]= C_inv * pt_homogeneous;
    }


}
void PnPSolver::ComputeM(const std::vector<Eigen::Vector2d> &pts2d, const Eigen::Matrix3d &K,
                         const std::vector<Eigen::Vector4d> &hbcs, Eigen::MatrixXd &M) {

    const int num_of_pts=pts2d.size();
    M.resize(2*num_of_pts,12);
    double f_x=K(0,0),f_y=K(1,1),c_x=K(0,2),c_y=K(1,2);

    for(int i=0;i<num_of_pts;i++){
        Eigen::Vector4d hb=hbcs[i];
        Eigen::Vector2d u=pts2d[i];
        for(int j=0;j<4;j++){
            M(2*i,j*3)=hb(j)*f_x;
            M(2*i,j*3+1)=0;
            M(2*i,j*3+2)=hb(j)*(c_x-u (0));

            M(2*i+1,j*3)=0;
            M(2*i+1,j*3+1)=hb(j)*f_y;
            M(2*i+1,j*3+2)=hb(j)*(c_y-u(1));

        }


    }




}
void PnPSolver::ComputeL6x10(const Eigen::Matrix<double, 12, 4> &V, Eigen::Matrix<double, 6, 10> &L6x10) {
    const int idx0[6] = {0, 0, 0, 1, 1, 2};
    const int idx1[6] = {1, 2, 3, 2, 3, 3};
    for(int i=0;i<6;i++){
        const int idi=idx0[i]*3;
        const int idj=idx1[i]*3;
        // the first control point.
        const Eigen::Vector3d v1i = V.block ( idi, 0, 3, 1 );
        const Eigen::Vector3d v2i = V.block ( idi, 1, 3, 1 );
        const Eigen::Vector3d v3i = V.block ( idi, 2, 3, 1 );
        const Eigen::Vector3d v4i = V.block ( idi, 3, 3, 1 );

        // the second control point
        const Eigen::Vector3d v1j = V.block ( idj, 0, 3, 1 );
        const Eigen::Vector3d v2j = V.block ( idj, 1, 3, 1 );
        const Eigen::Vector3d v3j = V.block ( idj, 2, 3, 1 );
        const Eigen::Vector3d v4j = V.block ( idj, 3, 3, 1 );

        Eigen::Vector3d S1 = v1i - v1j;
        Eigen::Vector3d S2 = v2i - v2j;
        Eigen::Vector3d S3 = v3i - v3j;
        Eigen::Vector3d S4 = v4i - v4j;

        Eigen::Matrix<double, 1, 3> S1_T = S1.transpose();
        Eigen::Matrix<double, 1, 3> S2_T = S2.transpose();
        Eigen::Matrix<double, 1, 3> S3_T = S3.transpose();
        Eigen::Matrix<double, 1, 3> S4_T = S4.transpose();

        //[B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
        L6x10 ( i, 0 ) = S1_T * S1;
        L6x10( i, 1 ) = 2 * S1_T * S2;
        L6x10( i, 2 ) = S2_T * S2;
        L6x10 ( i, 3 ) = 2 * S1_T * S3;
        L6x10( i, 4 ) = 2 * S2_T * S3;
        L6x10( i, 5 ) = S3_T * S3;
        L6x10 ( i, 6 ) = 2 * S1_T * S4;
        L6x10 ( i, 7 ) = 2 * S2_T * S4;
        L6x10 ( i, 8 ) = 2 * S3_T * S4;
        L6x10( i, 9 ) = S4_T * S4;


    }


}

void PnPSolver::ComputeRho(const std::vector<Eigen::Vector3d> &control_points_world, Eigen::VectorXd &rho) {
    rho.resize(6);

    const int idx0[6] = {0, 0, 0, 1, 1, 2};
    const int idx1[6] = {1, 2, 3, 2, 3, 3};
    for(int i=0;i<6;i++){
        Eigen::Vector3d dist=control_points_world[idx0[i]] -control_points_world[idx1[i]];
        rho(i)=dist.dot(dist);
    }

}


void PnPSolver::ComputeBetaApproximate1(const Eigen::Matrix<double, 6, 10> &L6x10, const Eigen::VectorXd &rho,
                                        Eigen::Vector4d &beta) {
    // betas10        = [B11 B12 B22 B13 B23 B33 B14 B24 B34 B44]
    // betas_approx_1 = [B11 B12     B13         B14]

    Eigen::Matrix<double ,6,4> L6x4;

    L6x4.col(0)=L6x10.col(0);
    L6x4.col(1)=L6x10.col(1);
    L6x4.col(2)=L6x10.col(3);
    L6x4.col(3)=L6x10.col(6);

    Eigen::VectorXd b4=L6x4.colPivHouseholderQr().solve(rho);

    if(b4(0)<0)
    {
        beta(0)=sqrt(-b4(0));
        beta(1)=-b4(1)/beta(0);
        beta(2)=-b4(2)/beta(0);
        beta(3)=(-b4(3))/beta(0);

    }
    else {
        beta(0)=sqrt(b4(0));
        beta(1)=b4(1)/beta(0);
        beta(2)=b4(2)/beta(0);
        beta(3)=(b4(3))/beta(0);

    }





}
void PnPSolver::ComputeJacobianMatrix(const Eigen::Matrix<double, 6, 10> &L6x10, const Eigen::Vector4d &beta,
                                      Eigen::Matrix<double, 6, 4> & jacobian) {
    Eigen::Matrix<double,10,4 > jacobian_beta_bata;
    double b1=beta(0),b2=beta(1),b3=beta(2),b4=beta(3);
    jacobian_beta_bata<<2*b1,0,0,0,
                        b2,b1,0,0,
                        0,2*b2,0,0,
                        b3,0,b1,0,
                        0,b3,b2,0,
                        0,0,2*b3,0,
                        b4,0,0,b1,
                        0,b4,0,b2,
                        0,0,b4,b3,
                        0,0,0,2*b4;
    jacobian=L6x10*jacobian_beta_bata;


}
void PnPSolver::ComputeResidual(const Eigen::Matrix<double, 6, 10> &L6x10, const Eigen::Vector4d &beta,
                                const Eigen::VectorXd &rho, Eigen::VectorXd &residual) {
    Eigen::VectorXd beta_vetor(10);
    double b1=beta(0),b2=beta(1),b3=beta(2),b4=beta(3);
    beta_vetor<<b1*b1,
                b1*b2,
                b2*b2,
                b1*b3,
                b2*b3,
                b3*b3,
                b1*b4,
                b2*b4,
                b3*b4,
                b4*b4;
    residual= L6x10 * beta_vetor - rho;

}


void PnPSolver::DoGaussNewtonOptimization(const Eigen::Matrix<double, 6, 10> &L6x10, const Eigen::VectorXd &rho,
                                           Eigen::Vector4d& beta) {
    const int iterattion_times=10;
    for(int i=0;i<iterattion_times;i++){

        Eigen::Matrix<double ,6,4> jacobian;
        Eigen::VectorXd residual(6);
        ComputeJacobianMatrix(L6x10,beta,jacobian);
        ComputeResidual(L6x10,beta,rho,residual);
        Eigen::Matrix<double , 4,6> jacobian_trans=jacobian.transpose();
        Eigen::Matrix4d H= jacobian_trans*jacobian;
        Eigen::Vector4d delta_beta =H.colPivHouseholderQr().solve(-jacobian_trans*residual);
        beta+=delta_beta;

    }



}
void PnPSolver::ComputeControlPointsCamera(const Eigen::Vector4d &beta, const Eigen::Matrix<double, 12, 4> &V,
                                           std::vector<Eigen::Vector3d>& control_points_camera) {
    control_points_camera.resize(4);
    Eigen::VectorXd  control_points_vector(12);
    control_points_vector=beta(0)*V.col(0)+beta(1)*V.col(1)+beta(2)*V.col(2)+beta(3)*V.col(3);

    for(int i=0;i<4;i++){
        control_points_camera[i]=control_points_vector.block(i*3,0,3,1); //control_points_vector.segment(i*3,3);

    }




}
void PnPSolver::ComputePts3DCameraCoordinates(const std::vector<Eigen::Vector4d>& hb_coordinates,
                                              const std::vector<Eigen::Vector3d> &control_points_camera,
                                              std::vector<Eigen::Vector3d>& pts_3d_camera) {
    const int num_of_pts=hb_coordinates.size();
    pts_3d_camera.resize(num_of_pts);
    for(int i=0;i<num_of_pts;i++){
        double a0=hb_coordinates[i](0),a1=hb_coordinates[i](1),a2=hb_coordinates[i](2),a3=hb_coordinates[i](3);
        pts_3d_camera[i]=a0*control_points_camera[0]+a1*control_points_camera[1]+a2*control_points_camera[2]+a3*control_points_camera[3];
    }

    for ( auto &pt:pts_3d_camera){
        if(pt(2)<0.0)
        {
            pt=-pt;
//        std::cout<<"pts_3d_camera z<0"<<std::endl;
//        std::cout<<pt<<std::endl;

        }
    }



}
void PnPSolver::ComputeRt(const std::vector<Eigen::Vector3d> &pts3d, const std::vector<Eigen::Vector3d> &pts_3d_camera,
                          Eigen::Matrix3d &R, Eigen::Vector3d& t) {
    const int num_of_pts=pts3d.size();

    Eigen::Vector3d sum_pts_world(0,0,0);

    for(const auto & pt:pts3d){
        sum_pts_world+=pt;
    }
    Eigen::Vector3d p_center_world=sum_pts_world/num_of_pts;


    Eigen::MatrixXd A(num_of_pts,3);
    for(int i=0;i<num_of_pts;i++){
        A.row(i)=(pts3d[i]-p_center_world).transpose();
    }

    Eigen::Vector3d sum_pts_camera(0,0,0);

    for(const auto & pt:pts_3d_camera){
        sum_pts_camera+=pt;
    }
    Eigen::Vector3d p_center_camera=sum_pts_camera/num_of_pts;

    Eigen::MatrixXd B(num_of_pts,3);
    for(int i=0;i<num_of_pts;i++){
        B.row(i)=(pts_3d_camera[i] - p_center_camera).transpose();

    }

    Eigen::Matrix3d H=B.transpose()*A;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();

    R = U*V.transpose();
    double detR = R.determinant();

    if (detR < 0){
     R.row(2)=-R.row(2);
    }

    t =p_center_camera  - R*p_center_world;



}