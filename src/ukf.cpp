#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;
  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  //TODO: when process noise sd are both 0.1 it's slow

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // TODO - do I need to initialize this here at all?

  // initial state vector
  x_ = VectorXd(n_x_);
  x_ << 1,1,1,1,1;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ <<   1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0 ,0, 1;

  //create augmented mean and covariance
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



  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  - Number of sigma points needs to be based on augmented number
  - Transition matrix F
  - Obs model - laser H_laser_
  - Obs model - radar H_radar_
  - Predicted covariance matrices for laser and radar R_laser and R_radar
  - Augmented x and P
  - predicted sigma points matrix
  - weights of sigma points

  TODO:
  One or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} measurement_pack The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.

  Determine detla_t and predict
  Update with either lidar or radar
  */

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
}

void UKF::InitializeFirstMeasurement(MeasurementPackage measurement_pack) {
  cout << "initialization" << endl;
  previous_timestamp_ = measurement_pack.timestamp_;
  if(measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //convert radar to CTRV to set first state
    cout << "INITIALIZE RADAR" << endl;
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
    cout << "INITIALIZE LASER" << endl;
    x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
  }
  x_aug_.head(5) = x_;
  is_initialized_ = true;
  //TODO: not clear when I need to update x_ vs x_aug_
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  GenerateSigmaPoints();
  PredictSigmaPoints(dt);
  PredictMeanCovariance();

  // TODO: figure out where X and X aug get updated.

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
  //weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  weights(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }

  //mean
  VectorXd x = VectorXd(n_x_);
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x = x+ weights(i) * Xsig_pred_.col(i);
  }

  //covariance
  MatrixXd P = MatrixXd(n_x_, n_x_);
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights(i) * x_diff * x_diff.transpose();
  }
  x_ = x;
  P_ = P;
  /*Update x, p
  TODO: do I need to update augmented x and P?
  Probably, b/c in real life you might not get measurements each time
  and need to predict further out, creating new sigma points?
  I seem to get different answers if I do or don't comment these out...
  */
  x_aug_.head(5) = x_;
  P_aug_.topLeftCorner(5,5) = P_;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} measurement_pack
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} measurement_pack
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

