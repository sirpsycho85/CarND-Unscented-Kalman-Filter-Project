#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/* TODO later:
clean up .h, public/private
figure out if/why the process noises that converged make physical sense
figure out why it's sensitive to initial x and P
*/

UKF::UKF() {

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  is_initialized_ = false;
  previous_timestamp_ = 0;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // initialize x, P
  x_ = VectorXd(n_x_);
  //x_ << 1,1,1,1,1;

  P_ = MatrixXd(n_x_, n_x_);
  P_ <<   1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0 ,0, 1;

  // augmented x, P
  x_aug_ = VectorXd(7);
  P_aug_ = MatrixXd(7, 7);

  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //create sigma point matrices
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  n_z_radar_ = 3;
  n_z_lidar_ = 2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  R_lidar_ = MatrixXd(n_z_lidar_,n_z_lidar_);
  R_lidar_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  R_radar_ = MatrixXd(n_z_radar_,n_z_radar_);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;
}


UKF::~UKF() {}


void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  if(!is_initialized_) {
    InitializeFirstMeasurement(measurement_pack);
  }

  double dt;
  dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  UKF::Prediction(dt);
  
  if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(measurement_pack);
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(measurement_pack);
  }
}


void UKF::InitializeFirstMeasurement(MeasurementPackage measurement_pack) {
  previous_timestamp_ = measurement_pack.timestamp_;
  if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //convert radar to CTRV to set first state
    float ro;
    float theta;
    float ro_dot;

    ro = measurement_pack.raw_measurements_[0];
    theta = measurement_pack.raw_measurements_[1];
    ro_dot = measurement_pack.raw_measurements_[2];

    float px;
    float py;
    float v;
    float psi;
    float psi_dot;

    px = ro * cos(theta);
    py = ro * sin(theta);
    v = ro_dot;
    psi = 0;
    psi_dot = 0;

    x_ << px, py, v, psi, psi_dot;
  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
  }
  x_aug_.head(5) = x_;
  is_initialized_ = true;
}


void UKF::Prediction(double dt) {
  GenerateSigmaPoints();
  PredictSigmaPoints(dt);
  PredictMeanCovariance();
}


void UKF::GenerateSigmaPoints() {
  A_ = P_aug_.llt().matrixL();
  
  Xsig_aug_.col(0)  = x_aug_;

  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_+n_aug_) * A_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * A_.col(i);
  }
}


void UKF::PredictSigmaPoints(double dt) {
  for (int i = 0; i< 2*n_aug_+1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*dt) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*dt) );
    }
    else {
        px_p = p_x + v*dt*cos(yaw);
        py_p = p_y + v*dt*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*dt;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*dt*dt * cos(yaw);
    py_p = py_p + 0.5*nu_a*dt*dt * sin(yaw);
    v_p = v_p + nu_a*dt;

    yaw_p = yaw_p + 0.5*nu_yawdd*dt*dt;
    yawd_p = yawd_p + nu_yawdd*dt;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}


void UKF::PredictMeanCovariance() {
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

  x_aug_.head(5) = x_;
  P_aug_.topLeftCorner(5,5) = P_;
  
  /* TODO: Pondering...do I need to update x_aug_ and P_aug_?
  In this project it has no impact.
  But in real life you might not get measurement for a while
  And need to do multiple predict steps in a row
  So you want the most up to date x_aug_ and P_aug_
  */
}


void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  //MEASUREMENT PREDICTION
  int n_z = n_z_lidar_;
  MatrixXd R = R_lidar_;
  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Zsig(0,i) = Xsig_pred_(0,i);  //px
    Zsig(1,i) = Xsig_pred_(1,i);  //py
  }

  UpdateCommon(n_z, R, Zsig, measurement_pack);
}


void UKF::UpdateRadar(MeasurementPackage measurement_pack) {

  //MEASUREMENT PREDICTION
  int n_z = n_z_radar_;
  MatrixXd R = R_radar_;
  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);              //r
    Zsig(1,i) = atan2(p_y,p_x);                       //phi
    if(sqrt(p_x*p_x + p_y*p_y) < 0.0001) {            //r_dot check for div by 0
      Zsig(2,i) = (p_x*v1 + p_y*v2) / 0.0001;
    }
    else {
      Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
    }
    
  }

  UpdateCommon(n_z, R, Zsig, measurement_pack);
}


void UKF::UpdateCommon(int n_z, MatrixXd R, MatrixXd Zsig, MeasurementPackage measurement_pack) {
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R;

  //MEASUREMENT UPDATE

  VectorXd z = measurement_pack.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state and augmented state
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  x_aug_.head(5) = x_;
  P_aug_.topLeftCorner(5,5) = P_;
}