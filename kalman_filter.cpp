#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
     TODO:
     * predict the state
     */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     TODO:
     * update the state by using Kalman Filter equations
     */
    
    /** MatrixXd H_KF = MatrixXd(2,4);
     H_KF = H_.topRows(2);
     VectorXd z_pred = H_KF * x_;
     VectorXd y = z - z_pred;
     MatrixXd Ht_KF = H_KF.transpose();
     MatrixXd S = H_KF * P_ * Ht_KF + R_.bottomRightCorner(2,2);
     MatrixXd Si = S.inverse();
     MatrixXd PHt = P_ * Ht_KF;
     MatrixXd K = PHt * Si;
     
     //new estimate
     x_ = x_ + (K * y);
     long x_size = x_.size();
     MatrixXd I = MatrixXd::Identity(x_size, x_size);
     P_ = (I - K * H_KF) * P_;*/
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
     TODO:
     * update the state by using Extended Kalman Filter equations
     */
    
    /**MatrixXd H_EKF =  MatrixXd(3,4);
     MatrixXd Hj = Tools::CalculateJacobian(const VectorXd x_); // used for calculating S,K,P
     H_.bottomRows(3) = Hj;
     H_EKF = H_.bottomRows(3);
     
     VectorXd h_x = Tools::CartisianToPolar(const VectorXd x_);
     VectorXd y = z - h_x;
     // adjusting the angle of polar residual measurement
     if(y[1] > (M_PI)){
     y[1] << y[1] - 2*M_PI;
     }
     else if(y[1] < (-M_PI)){
     y[1] << y[1] + 2*M_PI;
     }
     MatrixXd H_EKF_t = H_EKF.transpose();
     MatrixXd S = H_EKF * P_ * H_EKF_t + R_.topLeftCorner(3,3);
     MatrixXd Si = S.inverse();
     MatrixXd PHjt = P_ * H_EKF_t;
     MatrixXd K = PHjt * Si;
     
     //new estimate
     x_ = x_ + (K * y);
     long x_size = x_.size();
     MatrixXd I = MatrixXd::Identity(x_size, x_size);
     P_ = (I - K * H_EKF) * P_;*/
    
    
    // Cartisian to Polar
    VectorXd h_x(3);
    h_x <<  0.0, 0.0, 0.0;
    //recover state parameters
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);
    
    //pre-compute a set of terms to avoid repeated calculation
    double c1 = sqrt(px*px+py*py);
    double c2 = atan2(py,px);
    double c3 = (px*vx+py*vy)/c1;
    
    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "CartisianToPolar () - Error - Division by Zero" << endl;
    }
    
    //compute the Jacobian matrix
    h_x << c1, c2, c3;
    
    VectorXd y = z - h_x;
    // adjusting the angle of polar residual measurement
    if(y[1] > (M_PI)){
        y[1] = y[1] - 2*M_PI;
    }
    else if(y[1] < (-M_PI)){
        y[1] = y[1] + 2*M_PI;
    }
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
