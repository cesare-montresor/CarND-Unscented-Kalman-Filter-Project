#include "ukf.h"
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
  nis_lidar_ = ofstream("nis_lidar.csv");
  nis_radar_ = ofstream("nis_radar.csv");
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  is_initialized_ = false;
 

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
  
  
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;
  
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.5;
  
  
  ///* State dimension
  n_x_ = 5;
  
  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;
  
  sigma_count_ = 2 * n_aug_ + 1;
  
  
  
  sigma_spread_ = sqrt(lambda_+n_aug_);///* Sigma point spreading parameter
  
  ///* Weights of sigma points
  weights_ = VectorXd(sigma_count_);
  
  
  
  
  
  
  //init mean state vector
  x_ = VectorXd(n_x_);
  x_aug_ = VectorXd(n_aug_);
  
  //init covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_aug_ = MatrixXd::Identity(n_aug_, n_aug_);
  P_.diagonal()<<0,0,1,1,1;
  
  //init augmented mean state vector
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  
  //init augmented covariance matrix
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  //init weights
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<sigma_count_; i++) {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }
  
  R_radar_ = MatrixXd(3,3);
  R_radar_ <<
  std_radr_*std_radr_, 0, 0,
  0, std_radphi_*std_radphi_, 0,
  0, 0, std_radrd_*std_radrd_;
  
  
  R_lidar_ = MatrixXd(2,2);
  R_lidar_ <<
  std_laspx_*std_laspx_, 0,
  0, std_laspy_*std_laspy_;
  
}

UKF::~UKF() {
  nis_lidar_.close();
  nis_radar_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  
  if (!is_initialized_){
    VectorXd meas = meas_package.raw_measurements_;
    previous_timestamp_ = meas_package.timestamp_;
    if ( meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas[0],meas[1],0,0,0;
    }else if( use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR){
      x_<< t.Polar2Cartesian(meas[0], meas[1], meas[2]),0;
    }
    
    /*
    cout<<"\n---------------- AFTER INIT ----------------\n";
    cout<<"X:\n"<<x_<<"\n";
    cout<<"P:\n"<<P_<<"\n";
*/
    
    is_initialized_ = true;
    return;
  }
  
  
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  
  
  Prediction(dt);
  /*
  cout<<"\n---------------- AFTER PREDICT ----------------\n";
  cout<<"X:\n"<<x_<<"\n";
  cout<<"P:\n"<<P_<<"\n";
  */
  if ( use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar( meas_package );
  }else if( use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR){
    UpdateRadar( meas_package );
  }
  
  x_aug_.head(5) = x_;
  P_aug_.topLeftCorner(5,5) = P_;
  /*
  cout<<"\n---------------- AFTER UPDATE ----------------\n";
  cout<<"X:\n"<<x_<<"\n";
  cout<<"P:\n"<<P_<<"\n";
  */
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
/**
 TODO:
 
 Complete this function! Estimate the object's location. Modify the state
 vector, x_. Predict sigma points, the state, and the state covariance matrix.
 */


void UKF::Prediction(double delta_t) {

  MatrixXd Xsig_aug = GenerateSigmaPoints();
  Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t);
  UpdateMeanAndCovariance();
}


MatrixXd UKF::GenerateSigmaPoints(){
  
  ///* predicted sigma points matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, sigma_count_);
  
  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug_ + sigma_spread_ * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug_ - sigma_spread_ * L.col(i);
  }
  
  return Xsig_aug;
}

MatrixXd UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t){
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd::Zero(n_x_, sigma_count_);
  //predict sigma points
  for (int i = 0; i< sigma_count_; i++)
  {
    //extract values for better readability
    double p_x  = Xsig_aug(0,i);
    double p_y  = Xsig_aug(1,i);
    double v    = Xsig_aug(2,i);
    double yaw  = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
    
    //predicted state values
    double px_p, py_p;
    
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }
    
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    
    //add process noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    
    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
    //cout<<"Xsig_pred:\n"<<Xsig_pred.col(i)<<"\n\n";
  }
  
  return Xsig_pred;
}

void UKF::UpdateMeanAndCovariance(){
  x_.fill(0.0); // reset X
  P_.fill(0.0); // Reset P
  
  //predicted state mean
  for (int i = 0; i < sigma_count_; i++) {
    x_ = x_ + ( weights_(i) * Xsig_pred_.col(i) ); // Generate new means
  }
  
  //predicted state covariance matrix
  for (int i = 0; i < sigma_count_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_; // state difference
    x_diff(3) = t.NormalizeAngle( x_diff(3) );
    P_ = P_ + ( weights_(i) * x_diff * x_diff.transpose() );
  }
}




/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
  int n_z = 2;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, sigma_count_);
  
  //transform sigma points into measurement space
  for (int i = 0; i < sigma_count_ ; i++) {
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    
    // measurement model
    Zsig(0,i) = p_x;                        //x
    Zsig(1,i) = p_y;                        //y
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i=0; i < sigma_count_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i = 0; i < sigma_count_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;        //residual
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  S = S + R_lidar_; //add measurement noise covariance matrix
  
  
  
  //calculate cross correlation matrix
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < sigma_count_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;  //residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_; // state difference
    
    z_diff(1) = t.NormalizeAngle(z_diff(1));
    x_diff(3) = t.NormalizeAngle(x_diff(3));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  
  MatrixXd S_inverse = S.inverse();
  
  MatrixXd K = Tc * S.inverse(); //Kalman gain K;
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred; //residual
  z_diff(1) = t.NormalizeAngle(z_diff(1));
  
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  float nis = z_diff.transpose() * S_inverse * z_diff;
  nis_lidar_ << nis << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, sigma_count_);
  //transform sigma points into measurement space
  for (int i = 0; i < sigma_count_ ; i++) {  //2n+1 simga points
    cout<<Xsig_pred_.col(i);
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    double r = sqrt(p_x*p_x + p_y*p_y);
    
    // measurement model
    Zsig(0,i) = r;                                              //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    if (r > 0.001) { //avoid division by zero near the origin
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / r;                       //r_dot
    }
    else{
      //need a better way to initialize R_dot, it skerw my NIS and RMSE for the first 1-2 readings :(
      Zsig(2,i) = 0;
    }
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd::Zero(n_z);
  for (int i=0; i < sigma_count_; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i = 0; i < sigma_count_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;     //residual
    z_diff(1) = t.NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  //add measurement noise covariance matrix
  S = S + R_radar_;
  
  
  //calculate cross correlation matrix
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i = 0; i < sigma_count_; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;   //residual
    VectorXd x_diff = Xsig_pred_.col(i) - x_; // state difference
    
    z_diff(1) = t.NormalizeAngle(z_diff(1));
    x_diff(3) = t.NormalizeAngle(x_diff(3));
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  
  MatrixXd S_inverse = S.inverse();

  MatrixXd K = Tc * S_inverse; //Kalman gain K;
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred; //residual
  z_diff(1) = t.NormalizeAngle(z_diff(1));
//cout<<"\nraw_measurements_:\n"<<meas_package.raw_measurements_<<"\n";
//  cout<<"\nz_diff:\n"<<z_diff<<"\n";
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  // calculate NIS
  float nis = z_diff.transpose() * S_inverse * z_diff;
  nis_radar_ << nis << endl;
}

