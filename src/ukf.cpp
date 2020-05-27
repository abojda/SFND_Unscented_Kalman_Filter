#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */


  // Set state and augmented state dimensions
  n_x_ = 5;
  n_aug_= 7;

  time_us_ = 0;

  // Set sigma point spreading parameter
  lambda_= 3 - n_aug_;
  
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initialize weights_ vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  
  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
  {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // Create radar measurement noise matrix
  R_radar_ = MatrixXd(3, 3);
  R_radar_.fill(0);
  R_radar_(0, 0) = std_radr_ * std_radr_;
  R_radar_(1, 1) = std_radphi_ * std_radphi_;
  R_radar_(2, 2) = std_radrd_ * std_radrd_;

  // Create radar measurement noise matrix
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_.fill(0);
  R_lidar_(0, 0) = std_laspx_ * std_laspx_;
  R_lidar_(1, 1) = std_laspy_ * std_laspy_;

  // Initialize state covariance
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_(0,0) = std_laspx_*std_laspy_;
  P_(1,1) = std_laspx_*std_laspy_;

  // NIS printing
  print_NIS_lidar_ = false;
  print_NIS_radar_ = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      double x = meas_package.raw_measurements_[0];
      double y = meas_package.raw_measurements_[1];
      x_ << x, y, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      x_ << x, y, 0, 0, 0;
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // Prediction step
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  
  // Measurement step
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  VectorXd x_aug_ = VectorXd(n_aug_);
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  P_aug_.fill(0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  // Create square root matrix
  MatrixXd A = P_aug_.llt().matrixL();
  A *= sqrt(lambda_ + n_aug_);
  
  // Generate sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; ++i)
  {
      Xsig_aug_.col(i + 1) = x_aug_ + A.col(i);
      Xsig_aug_.col(i + n_aug_ + 1) = x_aug_ - A.col(i);
  }
  
  // Predict sigma points
  double x, y, v, psi, psi_d, nu_a, nu_psi_dd;
  double x_k1, y_k1, v_k1, psi_k1, psi_d_k1;

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    x = Xsig_aug_.col(i)[0];
    y = Xsig_aug_.col(i)[1];
    v = Xsig_aug_.col(i)[2];
    psi = Xsig_aug_.col(i)[3];
    psi_d = Xsig_aug_.col(i)[4];
    nu_a = Xsig_aug_.col(i)[5];
    nu_psi_dd = Xsig_aug_.col(i)[6];
    
    // Deterministic part
    if (fabs(psi_d) > 0.001)
    {
        // Normal case
        x_k1 = x + v / psi_d * (sin(psi + psi_d * delta_t) - sin(psi));
        y_k1 = y + v / psi_d * (-cos(psi + psi_d * delta_t) + cos(psi));
    }
    else
    {
        // Special case to avoid division by zero
        x_k1 = x + v * cos(psi)*delta_t;
        y_k1 = y + v * sin(psi)*delta_t;
    }

    v_k1 = v + 0;
    psi_k1 = psi + (psi_d * delta_t);
    psi_d_k1 = psi_d + 0;

    // Stochastic (noise) part
    x_k1 += cos(psi) * nu_a * delta_t*delta_t / 2.0;
    y_k1 += sin(psi) * nu_a * delta_t*delta_t / 2.0;
    v_k1 += nu_a * delta_t;
    psi_k1 += nu_psi_dd * delta_t*delta_t / 2.0;
    psi_d_k1 += nu_psi_dd * delta_t;

    // Put into the correct column
    Xsig_pred_.col(i) << x_k1, y_k1, v_k1, psi_k1, psi_d_k1;
  }
 
  VectorXd new_x_ = VectorXd(n_x_);
  MatrixXd new_P_ = MatrixXd(n_x_, n_x_);

  new_x_.fill(0);
  new_P_.fill(0);

  // Calculate state mean
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      new_x_ += (weights_[i] * Xsig_pred_.col(i));
  }

  // Calculate state covariance matrix
  for (int i = 0; i < 2*n_aug_+1; ++i)
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    new_P_ += weights_[i] * x_diff * x_diff.transpose();
  }

  x_ = new_x_;
  P_ = new_P_;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  // Measurement dimension
  int n_z_ = 2;

  VectorXd Z = meas_package.raw_measurements_;
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd Z_pred = VectorXd(n_z_);
  MatrixXd S_ = MatrixXd(n_z_, n_z_);

  Z_pred.fill(0);
  S_.fill(0);

  // Transform sigma points into measurement space - reuse previous sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      double x = Xsig_pred_(0, i);
      double y = Xsig_pred_(1, i);
      Zsig_.col(i) << x, y; 
  }
  
  // calculate mean predicted measurement
  for (int i = 0; i< 2 * n_aug_ + 1; ++i)
  {
    Z_pred += (weights_[i] * Zsig_.col(i));
  }
  
  // calculate innovation covariance matrix S
  VectorXd z_diff = VectorXd(3);
  for (int i = 0; i< 2 * n_aug_ + 1; ++i)
  {
    z_diff = (Zsig_.col(i) - Z_pred);
    S_ += weights_[i] * z_diff * z_diff.transpose();
  }

  S_ += R_lidar_;

  // UKF Update
  VectorXd x_diff = VectorXd(n_x_);
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    x_diff = Xsig_pred_.col(i) - x_;
    z_diff = Zsig_.col(i) - Z_pred;

    Tc += weights_[i] * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  // update state mean and covariance matrix
  z_diff = Z - Z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ += K * z_diff;
  P_ -= K * S_ * K.transpose();

  // Calculate NIS for lidar
  NIS_lidar_ = z_diff.transpose() * S_.inverse() * z_diff;

  if (print_NIS_lidar_)
  {
    std::cout << NIS_lidar_ << std::endl;
  }
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

  // Measurement dimension
  int n_z_ = 3;

  VectorXd Z = meas_package.raw_measurements_;
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  VectorXd Z_pred = VectorXd(n_z_);
  MatrixXd S_ = MatrixXd(n_z_, n_z_);

  Z_pred.fill(0);
  S_.fill(0);

  // Transform sigma points into measurement space - reuse previous sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
      double x = Xsig_pred_(0, i);
      double y = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);
      double psi_d = Xsig_pred_(4, i);
      
      double rho = sqrt(x*x + y*y);
      double phi = atan2(y, x);
      // double rho_d = (x * cos(psi) * v + y * sin(psi) * v) / sqrt(x*x + y*y);
      double rho_d = (x * cos(psi) * v + y * sin(psi) * v) / std::max(0.0001, rho);
      

      Zsig_.col(i) << rho, phi, rho_d; 
  }
  
  // calculate mean predicted measurement
  for (int i = 0; i< 2 * n_aug_ + 1; ++i)
  {
    Z_pred += (weights_[i] * Zsig_.col(i));
  }
  
  // calculate innovation covariance matrix S
  VectorXd z_diff = VectorXd(3);
  for (int i = 0; i< 2 * n_aug_ + 1; ++i)
  {
    z_diff = (Zsig_.col(i) - Z_pred);

    // Angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S_ += weights_[i] * z_diff * z_diff.transpose();
  }

  S_ += R_radar_;


  // UKF Update
  VectorXd x_diff = VectorXd(n_x_);
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0);

  // calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    x_diff = Xsig_pred_.col(i) - x_;
    z_diff = Zsig_.col(i) - Z_pred;

    // Angle normalization    
    while (x_diff(3) > M_PI)    x_diff(3) -= 2.*M_PI;
    while (x_diff(3) <-M_PI)    x_diff(3) += 2.*M_PI;
    while (z_diff(1) > M_PI)    z_diff(1) -= 2.*M_PI;
    while (z_diff(1) <-M_PI)    z_diff(1) += 2.*M_PI;

    Tc += weights_[i] * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  // update state mean and covariance matrix
  z_diff = Z - Z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ += K * z_diff;
  P_ -= K * S_ * K.transpose();

  // Calculate NIS for radar
  NIS_radar_ = z_diff.transpose() * S_.inverse() * z_diff;

  if (print_NIS_radar_)
  {
    std::cout << NIS_radar_ << std::endl;
  }
}
