#ifndef ESTIMATOR_H
#define ESTIMATOR_H

namespace bifoiler {

class MEKF<typename Dynamics, typename Observation> {
public:
    using Scalar = typename Dynamics::State::Scalar;
    using State = Eigen::Matrix<Scalar, 18, 1>; // [v, w, r, a, bg, ba];
    using SysState = Eigen::Matrix<Scalar, 13, 1>; // [v, w, r, q_BI];
    using Control = typename Dynamics::Control;
    using Measurement = typename Observation::Measurement;
    using StateCov = Eigen::Matrix<Scalar, State::RowsAtCompileTime, State::RowsAtCompileTime>;
    using Quaternion = Eigen::Matrix<Scalar, 4, 1>; // [w, x, y, z];
    using Vector3 = const Eigen::Matrix<Scalar, 3, 1>;
    enum {
        nx = State::RowsAtCompileTime,
        nz = Measurement::RowsAtCompileTime,
    };

private:
    Dynamics f;
    Observation h;
    State x; // state vector
    Quaternion qref;
    StateCov P; // state covariance
    Eigen::Matrix<Scalar, nz, nz> S; // innovation covariance
    Eigen::Matrix<Scalar, nx, nz> K; // Kalman gain
    Eigen::Matrix<Scalar, nx, nx> I; // identity

    Quaternion _qerr(const Vector3 &a);

public:
    SysState get_system_state();

    void predict(const Control &u);
    void correct(const Measurement &z);
    void update(const Control &u, const Measurement &z);

};

} // namespace bifoiler

#endif /* ESTIMATOR_H */