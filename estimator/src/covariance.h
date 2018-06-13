// from get_cov.m "State Estimation for a Hydrofoil Boat" by Martin Schwendener

#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <Eigen/Dense>
#include <cmath>

#define DEG2RAD(x) ((M_PI / 180.0f) * (x))

// Guesses of Process noise and Measurement Noise for the MTi-G-710 by XSENS

template <typename Scalar>
static Eigen::Matrix<Scalar, 18, 18> process_noise_covariance(void)
{
    const Scalar t_samp = 0.02; // TODO: parametrize

    Eigen::Matrix<Scalar, 18, 18> Q;
    Eigen::Matrix<Scalar, 3, 1> Qv, Qw,Qr, Qa, Qbg, Qba;

    Q.setZero();

    // Process Noise/Disturbances
    // ---> How far  can i be of in 0.02 seconds? (2 cm/0.02s = 1m/s)

    // Allthough the process noise appears in the RHS, it has the same unit as
    // the state, since Difference Equations!

    // Velocity Disturbance. unaccounted kinetics
    Qv << 0.01, 0.8, 0.4;             // [(m/s)^2], Wild guess

    // Angular Vel Disturbance. unaccounted kinetics
    //Qw << 0.025, 0.025, 0.1;          // (rad/s)^2, Wild guess
    Qw << 0.1, 0.8, 1.2;              // [(rad/s)^2]


    // Position Disturbance. Emerge by the integration of kinetics that we don't model
    Qr << 0.05, 0.05, 0.03;           // [m^2], Wild guess

    // Attitude Disturbance
    Qa << 0.005, 0.005, 0.005;        // [rad^2], Wild guess

    // Gyro bias stability/variance (Change of Bias with time, process noise of
    // noise state)
    Scalar bs_gyro = 10;                        // [deg/h]
    Scalar sigma_b_gyro = DEG2RAD(bs_gyro)/3600;// [rad/s]
    Qbg.setOnes();
    Qbg *= pow(sigma_b_gyro, 2);

    // Accelerometer Bias stability/variance
    Scalar bs_acc = 40;                             // [mug]
    Scalar sigma_b_acc = bs_acc/1e6*9.81*t_samp;    // [m/s]
    Qba.setOnes();
    Qba *= pow(sigma_b_acc, 2);     // [(m/s)^2]

    Q.setZero();
    Q.diagonal() << Qv, Qw,Qr, Qa, Qbg, Qba;  // Process noise

    return Q;
}

template <typename Scalar>
static Eigen::Matrix<Scalar, 12, 12> measurement_noise_covariance(void)
{
    const Scalar t_samp = 0.02; // TODO: parametrize

    // Measurement Noise
    Eigen::Matrix<Scalar, 12, 12> R;
    Eigen::Matrix<Scalar, 3, 1> Rv, Rbg, Rr, Rba;

    // GNS Variance
    Rr << pow(2.5, 2), pow(2.5, 2), pow(5.0, 2); // [m^2]

    Rv.setOnes();
    Rv *= pow(0.05, 2);                    // [(m/s)^2], std dev @ 30m/s

    // Gyro variance, noise (not density, since discrete time equation)
    Scalar rho_sigma_gyro = 0.01; // [deg/s/sqrt(Hz)]
    Scalar sigma_gyro = DEG2RAD(rho_sigma_gyro) * sqrt(1.0f/t_samp); // [rad/s]
    Rbg.setOnes();
    Rbg *= pow(sigma_gyro, 2);

    // Accelerometer Variance
    // Noise density. Max = 150, typical: 80
    Scalar rho_sigma_acc = 150; // [mug/sqrt(Hz)]
    Scalar sigma_acc = rho_sigma_acc / 1e6*9.81*t_samp*sqrt(1.0f/t_samp); // [m/s]

    Rba.setOnes();
    Rba *= pow(sigma_acc, 2); // [(m/s)^2]

    R.setZero();
    R.diagonal() << Rv, Rbg, Rr, Rba;

    return R;
}

#endif /* COVARIANCE_H */
