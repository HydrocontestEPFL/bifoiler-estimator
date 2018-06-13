// from get_cov.m "State Estimation for a Hydrofoil Boat" by Martin Schwendener

#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <Eigen/Dense>
#include <math.h>

// Guesses of Process noise and Measurement Noise for the MTi-G-710 by XSENS

template <typename Scalar>
static Eigen::Matrix<Scalar, 18, 18> process_noise_covariance(void)
{
    Eigen::Matrix<Scalar, 18, 18> Q;

    // Process Noise/Disturbances
    // ---> How far  can i be of in 0.02 seconds? (2 cm/0.02s = 1m/s)

    // Allthough the process noise appears in the RHS, it has the same unit as
    // the state, since Difference Equations!
#if 0
    // Velocity Disturbance. unaccounted kinetics
    Qv          = diag([0.01 0.8 0.4]);           // (m/s)^2, Wild guess

    // Angular Vel Disturbance. unaccounted kinetics
    //Qw          = diag([0.025 0.025 0.1]);     // (rad/s)^2, Wild guess
    Qw          = diag([0.1 0.8 1.2]);


    // Position Disturbance. Emerge by the integration of kinetics that we don't model
    Qr          = diag([0.05 0.05 0.03]);         // (m)^2, Wild guess

    // Attitude Disturbance
    Qa          = diag([0.005 0.005 0.005]);          // (rad)^2, Wild guess

    // Gyro bias stability/variance (Change of Bias with time, process noise of
    // noise state)
    bs          = 10;                               // deg/h
    sigma_b     = deg2rad(bs)/3600;                 // rad/s
    Qbg         = sigma_b^2*eye(3);

    // Accelerometer Bias stability/variance
    bs          = 40;                               // mug
    sigma_b     = bs/1e6*9.81*t_samp;               // m/s
    Qba         = sigma_b^2*eye(3);                 // m^2/s^2

    // Process Covariance Matrix
    Q           = blkdiag(Qv, Qw,Qr, Qa, Qbg, Qba);  // Process noise
#endif
    Q.setZero();
    Q.diagonal().setOnes(); // TODO: only for testing...

    return Q;
}

template <typename Scalar>
static Eigen::Matrix<Scalar, 12, 12> measurement_noise_covariance(void)
{
    // Measurement Noise
    Eigen::Matrix<Scalar, 12, 12> R;

#if 0
    // GNS Variance
    Rr          = diag([2.5^2 2.5^2 5^2]);          // m^2

    Rv          = 0.05^2*eye(3);                    // (m/s)^2, std dev @ 30m/s

    // Gyro variance, noise (not density, since discrete time equation)
    rho_sigma   = 0.01;                             // deg/s/sqrt(Hz)
    sigma_gyro  = deg2rad(rho_sigma)*sqrt(1/t_samp);// rad/s
    Rbg          = sigma_gyro^2*eye(3);

    // Accelerometer Variance
    // Noise density. Max = 150, typical: 80
    Scalar rho_sigma = 150;                              // mug/sqrt(Hz)
    Scalar sigma_acc = rho_sigma / 1e6*9.81*t_samp*sqrt(10.f/t_samp);     // m/s

    Eigen::Matrix<Scalar, 3, 1> Rba
    Rba.setOnes();
    Rba *= sigma_acc^2*eye(3); // m^2/s^2

    R.setZero();
    R.diagonal() << Rv, Rbg, Rr, Rba;
#endif

    R.setZero();
    R.diagonal().setOnes(); // TODO: only for testing...

    return R;
}

#endif /* COVARIANCE_H */
