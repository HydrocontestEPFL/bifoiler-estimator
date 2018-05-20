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
    xs.block<4,1>(9,0) = qref;
    xs.block<9,1>(0,0) = x.block<9,1>(0,0);
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
    SysState x_sys, x_sys_p
    x_sys = get_system_state();

    // numerical integration
    x_sys_p = Dynamics.integrate(x_sys, u);

    F = f.jacobian(x, u, qref);

    State xe_p;
    xe_p.block<9,1>(0,0) = x_sys_p.block<9,1>(0,0);
    xe_p.block<9,1>(9,0) = F.block<9,18>(9,0)*x;

    // Covariance prediction
    StateCov P_p = F * P * F.transpose() + f.Q;
}

void MEKF::correct(const Measurement &z)
{
    Measurement y;                      // innovation
    Eigen::Matrix<Scalar, nz, nx> H;    // jacobian of h
    Eigen::Matrix<Scalar, nx, nx> IKH;  // temporary matrix

    H = h.jacobian(xe_p, u, qref);
    S = H * P_p * H.transpose() + h.R;

    // efficiently compute: K = P * H.transpose() * S.inverse();
    K = S.llt().solve(H * P_p).transpose();

    y = z - h(x);
    IKH = (I - K * H);

    // Measurement update
    x = xe_p + K * y;
    P = IKH * P_p * IKH.transpose() + K * h.R * K.transpose();

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
