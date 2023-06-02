#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::pow;
double normalize_angle(double angle) {
	angle = std::fmod(angle + M_PI, 2 * M_PI);
	if (angle < 0)angle += 2 * M_PI;
	return angle - M_PI;
}

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
  std_a_ = 3;

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
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2 * n_aug_ + 1;
  weights_ = VectorXd(n_sig_);
  lambda_ = 3.0 - n_x_;
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  Xsig_pred_.fill(0.0);
  double w0 = lambda_ / (lambda_ + n_aug_);
  double w = 1 / (2 * (n_aug_ + lambda_));
  weights_.fill(w);
  weights_(0) = w0;
  is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::InitializeUKF(MeasurementPackage meas_package) {
	x_.fill(0.0);
	x_.head(2) << meas_package.raw_measurements_;
	time_us_ = meas_package.timestamp_;
	is_initialized_ = true;
}
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
	if (!is_initialized_) {
		InitializeUKF(meas_package);
		if (MeasurementPackage::SensorType::LASER == meas_package.sensor_type_) {
			UpdateLidar(meas_package);
		}
		else if (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_) {
			UpdateRadar(meas_package);
		}
	}
	else {
		double delta_t = (meas_package.timestamp_ - time_us_) / 1e6;
		time_us_ = meas_package.timestamp_;
		Prediction(delta_t);
		if (MeasurementPackage::SensorType::LASER == meas_package.sensor_type_) {
			UpdateLidar(meas_package);
		}
		else if (MeasurementPackage::SensorType::RADAR == meas_package.sensor_type_) {
			UpdateRadar(meas_package);
		}
	}
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
	//Generate sigma points
	MatrixXd Xsig = MatrixXd(n_aug_, n_sig_);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug(n_x_, n_x_) = pow(std_a_, 2);
	P_aug(n_x_ + 1, n_x_ + 1) = pow(std_yawdd_, 2);
	//Take square root of P matrix
	MatrixXd A = P_aug.llt().matrixL();
	double c_aug = sqrt(lambda_ + n_aug_);
	MatrixXd C_aug = c_aug*A;
	Xsig.col(0) = x_aug;
	for (int i = 1; i <= n_aug_; i++) {
		Xsig.col(i) = x_aug + C_aug.col(i - 1);
	}
	for (int i = n_aug_ + 1; i <= n_sig_ - 1; i++) {
		Xsig.col(i) = x_aug - C_aug.col(i - 1 - n_aug_);
	}
	//Predict sigma point
	double dt = delta_t;
	Xsig_pred_.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		double px = Xsig(0, i);
		double py = Xsig(1, i);
		double v = Xsig(2, i);
		double psi = Xsig(3, i);
		double psid = Xsig(4, i);
		double nu_a = Xsig(5, i);
		double nu_psidd = Xsig(6, i);

		// Transform from the augmented to the original (n_x-dim) state space

		// Avoid division by zero
		if (std::fabs(psid) > 0.001) {
			Xsig_pred_(0, i) = px + v / psid * (sin(psi + psid * dt) - sin(psi)) +
				pow(dt, 2) / 2 * cos(psi) * nu_a;
			Xsig_pred_(1, i) = py + v / psid * (-cos(psi + psid * dt) + cos(psi)) +
				pow(dt, 2) / 2 * sin(psi) * nu_a;
		}
		else {
			Xsig_pred_(0, i) = px + v * dt * cos(psi) +
				pow(dt, 2) / 2 * cos(psi) * nu_a;
			Xsig_pred_(1, i) = py + v * dt * sin(psi) +
				pow(dt, 2) / 2 * sin(psi) * nu_a;
		}

		Xsig_pred_(2, i) = v + 0 + dt * nu_a;
		Xsig_pred_(3, i) = psi + psid * dt + pow(dt, 2) / 2 * nu_psidd;
		Xsig_pred_(4, i) = psid + 0 + dt * nu_psidd;
	}
	VectorXd x = VectorXd(n_x_);
	x.fill(0.0);
	MatrixXd P = MatrixXd(n_x_, n_x_);
	P.fill(0.0);
	//Compute mean and covariance
	for (int i = 0; i < n_sig_; i++) {
		x = x + weights_(i)*Xsig_pred_.col(i);
	}
	for (int i = 0; i < n_sig_; i++) {
		p = P + weights_(i)*(Xsig_pred_.col(i) - x_)*(Xsig_pred_.col(i) - x_).transpose();
	}
	x_ = x;
	P_ = P;

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
	VectorXd z = meas_package.raw_measurements_;
	int n_z = z.size();
	MatrixXd H = MatrixXd(n_z, n_x_);
	H << 1, 0, 0, 0, 0,
		0, 1, 0, 0, 0;
	MatrixXd R = MatrixXd(n_z, n_z);
	R << pow(std_laspx_, 2), 0,
		0, pow(std_laspy_, 2);
	VectorXd z_pred = VectorXd(n_z);
	z_pred = x_.head(n_z);
	VectorXd y = z - z_pred;
	MatrixXd Ht = H.transpose();
	MatrixXd S = H*P_*Ht + R;
	MatrixXd Sinv = S.inverse();
	MatrixXd  Pht = P_*Ht;
	MatrixXd K = Pht*Sinv;
	MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
	x_ = x_ + (K*y);
	P_ = (I - K*H)*P_;

}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
	VectorXd z = meas_package.raw_measurements_;
	int n_z = z.size();
	MatrixXd Zsig = MatrixXd(n_z, n_sig_);
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	//transform the points to measurement space
	for (int i = 0; i < n_sig_; i++) {
		double px = Xsig_pred_(0, i);
		double py = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double psi = Xsig_pred_(3, i);
		double psid = Xsig_pred_(4, i);

		double rho = sqrt(pow(px, 2) + pow(py, 2));
		double phi = std::atan2(py, px);
		double rhod = 0.0;
		if (std::fabs(rho) > 0.001) {
			rhod = (px * cos(psi) * v + py * sin(psi) * v) / rho;
		}

		Zsig(0, i) = rho;
		Zsig(1, i) = phi;
		Zsig(2, i) = rhod;
	}
	//Compute mean and covariance of these points in measurement space
	for (int i = 0; i < n_sig_; i++) {
		z_pred = z_pred + weights_(i)*Zsig.col(i);
	}
	MatrixXd R = MatrixXd(n_z, n_z);
	double std_rho = pow(std_radr_, 2);
	double std_phi = pow(std_radphi_, 2);
	double std_rate = pow(std_radrd_, 2);
	double angle = 0.0;
	R << std_rho, 0, 0,
		0, std_phi, 0,
		0, 0, std_rate;
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		VectorXd zdiff = Zsig.col(i) - z_pred;
		zdiff(1) = normalize_angle(zdiff(1));
		S = S + weights_(i)*zdiff*zdiff.transpose();
	}
	S = S + R;
	//Compute Cross correlation
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);
	for (int i = 0; i < n_sig_; i++) {
		VectorXd zdiff = Zsig.col(i) - z_pred;
		zdiff(1) = normalize_angle(zdiff(1));
		VectorXd xdiff = Xsig_pred_.col(i) - x_;
		xdiff(3) = normalize_angle(xdiff(3));
		Tc += weights_(i)*xdiff*zdiff.transpose();
	}
	//Kalman gain
	MatrixXd K = MatrixXd(n_x_, n_z);
	MatrixXd Sinv = S.inverse();
	VectorXd y = z - z_pred;
	y(1) = normalize_angle(y(1));
	K = Tc*Sinv;
	x_ = x_ + K*y;
	MatrixXd Kt = K.transpose();
	P_ = P_ - K*S*Kt;


}