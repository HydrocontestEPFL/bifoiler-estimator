#include <Eigen/Dense>

namespace bifoiler {

class MEKF<typename Dynamics, typename Observation> {
public:
    using Scalar = typename Scalar;
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

    Quaternion _quat_error_mult(const Vector3 &a,
                                const Quaternion &qr);

public:
    SysState get_system_state();
    void update(const Control &u, const Measurement &z);

};

MEKF::SysState MEKF::get_system_state()
{
    SysState xs;
    xs << x.block<9,1>(0,0), qref;
    return xs;
}

MEKF::Quaternion MEKF::_quat_error_mul(const MEKF::Vector3 &a,
                                        const MEKF::Quaternion &qr)
{
    Quaternion dq;
    Quaternion q;

    // Factor
    const Scalar f = 1;

    // Scalar Part
    Scalar eta = 1 - a.squaredNorm() / 8;

    // Vector Part
    Vector3 nu = a / 2;

    // Error Parametrisation
    dq << eta, nu;
    dq *= f;

    // Erroneous Quaternion
    q = quatmul(dq, qr);

    return q;
}

void MEKF::MEKF()
{
    I.setIdentity();
}

void MEKF::predict(const Control &u)
{
    Eigen::Matrix<Scalar, nx, nx> F;
    SysState x_sys;

    x_sys = get_system_state();

    // numerical integration
    x_sys = Dynamics.integrate(x_sys, u);

    F = f.jacobian(x, u, qref);

    x << x_sys.block<9,1>(0,0), F.block<9,18>(9,0)*x

    // Covariance prediction
    P = F * P * F.transpose() + f.Q;
}

void MEKF::correct(const Measurement &z)
{
    Measurement y;                      // innovation
    Eigen::Matrix<Scalar, nz, nx> H;    // jacobian of h
    Eigen::Matrix<Scalar, nx, nx> IKH;  // temporary matrix

    H = h.jacobian(x, u, qref);
    S = H * P * H.transpose() + h.R;

    // efficiently compute: K = P * H.transpose() * S.inverse();
    K = S.llt().solve(H * P).transpose();

    y = z - h(x);
    IKH = (I - K * H);

    // Measurement update
    x = x + K * y;
    P = IKH * P * IKH.transpose() + K * h.R * K.transpose();

    // Attitude error transfer to reference quaternion qref
    Matrix<Scalar, 3, 1> a;
    a = x.block<3,1>(9,0);
    qref = _quat_error_mul(a, qref);

    // enforce unit norm constraint
    qref /= qref.norm();

    // reset MEKF error angle
    x.block<3,1>(9,0).setZero();
}

void MEKF::update(const Control &u, const Measurement &z)
{
    predict(u);
    correct(z);
}

} // namespace bifoiler