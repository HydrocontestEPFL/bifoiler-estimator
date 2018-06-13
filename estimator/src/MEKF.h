#ifndef ESTIMATOR_H
#define ESTIMATOR_H

namespace bifoiler {

class MEKF<typename Dynamics, typename Observation> {
public:
    using Scalar = typename Dynamics::EstimatorState::Scalar;
    using EstimatorState = Eigen::Matrix<Scalar, 18, 1>; // [v, w, r, a, bg, ba];
    using SystemState = Eigen::Matrix<Scalar, 13, 1>; // [v, w, r, q];
    using Control = typename Dynamics::Control;
    using Measurement = typename Observation::Measurement;
    using StateCov = Eigen::Matrix<Scalar, 18, 18>;
    using Quaternion = Eigen::Matrix<Scalar, 4, 1>; // [w, x, y, z];
    using Vector3 = const Eigen::Matrix<Scalar, 3, 1>;
    enum {
        nx = EstimatorState::RowsAtCompileTime,
        nxs = SystState::RowsAtCompileTime,
        nz = Measurement::RowsAtCompileTime,
    };

private:
    Dynamics f;
    Observation h;
    EstimatorState x; // state vector
    Quaternion qref;
    StateCov P; // state covariance
    Eigen::Matrix<Scalar, nz, nz> S; // innovation covariance
    Eigen::Matrix<Scalar, nx, nz> K; // Kalman gain
    Eigen::Matrix<Scalar, nx, nx> I; // identity

    Quaternion quatmul(const Quaternion &q1, const Quaternion &q2);
    Quaternion quat_error_mult(const Vector3 &a);

public:
    SystemState get_system_state();
    EstimatorState get_system_state() { return x; }
    StateCov get_state_covariance() { return P; }

    MEKF(const SystemState &x0, const StateCov &P0);

    void predict(const Control &u);
    void correct(const Measurement &z);
    void update(const Control &u, const Measurement &z);

};

} // namespace bifoiler

#endif /* ESTIMATOR_H */