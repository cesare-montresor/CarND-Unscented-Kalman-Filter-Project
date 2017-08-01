#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

#include "measurement_package.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  Tools t;
  ofstream nis_lidar_= ofstream("nis_lidar.csv");;
  ofstream nis_radar_= ofstream("nis_radar.csv");;

  
  // if this is false, laser measurements will be ignored (except during init)
  bool use_laser_ = true;
  
  // if this is false, radar measurements will be ignored (except during init)
  bool use_radar_ = true;
  
  bool is_initialized_ = false;
  
  
  // Laser measurement noise standard deviation position1 in m
  double std_laspx_ = 0.15;
  
  // Laser measurement noise standard deviation position2 in m
  double std_laspy_ = 0.15;
  
  // Radar measurement noise standard deviation radius in m
  double std_radr_ = 0.3;
  
  // Radar measurement noise standard deviation angle in rad
  double std_radphi_ = 0.03;
  
  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ = 0.3;
  
  
  
  int n_x_ = 5; ///* State dimension
  int n_aug_ = n_x_ + 2; ///* Augmented state dimension
  double lambda_ = 3 - n_aug_;
  double sigma_count_ = 2 * n_aug_ + 1;
  double sigma_spread_ = sqrt(lambda_+n_aug_); ///* Sigma point spreading parameter
  
  
  
  
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_;
  MatrixXd P_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;
  
  long long previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Weights of sigma points
  VectorXd weights_;
  
  

  MatrixXd R_radar_;
  MatrixXd R_lidar_;
  
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);
  
  // prediction sub operations
  MatrixXd GenerateSigmaPoints();
  MatrixXd PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t);
  void UpdateMeanAndCovariance();
  void ConvertSigmaToMesurment(int n_z, MatrixXd &S, MatrixXd &R);
  
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
};

#endif /* UKF_H */
