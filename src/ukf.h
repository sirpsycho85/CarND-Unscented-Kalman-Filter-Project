#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {

public:

  bool is_initialized_;
  long previous_timestamp_;

  // state vector [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  VectorXd x_aug_;

  // state covariance matrix
  MatrixXd P_;
  MatrixXd P_aug_;

  // square root of P_aug_
  MatrixXd A_;

  // sigma points matrices and weights
  MatrixXd Xsig_aug_;
  MatrixXd Xsig_pred_;
  VectorXd weights_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Measurement noise covariance matrices
  MatrixXd R_radar_;
  MatrixXd R_lidar_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  int n_z_radar_;
  int n_z_lidar_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  double NIS_radar_count_;
  double NIS_radar_sum_;
  double NIS_radar_mean_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();
  
  void InitializeFirstMeasurement(MeasurementPackage measurement_pack);
  void GenerateSigmaPoints();
  void PredictSigmaPoints(double dt);
  void PredictMeanCovariance();
  void PredictRadarMeasurement();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param dt Time between k and k+1 in s
   */
  void Prediction(double dt);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
  void UpdateCommon(int n_z, MatrixXd R, MatrixXd Zsig, MeasurementPackage measurement_pack);
};

#endif /* UKF_H */
