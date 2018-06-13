#include <iostream> // TODO: remove later

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "MEKF.h"

namespace bifoiler {

MEKF::SystemState MEKF::get_system_state()
{
    SystemState xs;
    xs << x.block<9,1>(0,0), qref;
    return xs;
}

MEKF::Quaternion MEKF::quatmul(const MEKF::Quaternion &q1, const MEKF::Quaternion &q2)
{
    Quaternion q;

    Scalar s1 = q1(0);
    Vector3 v1 = q1.block<3,1>(1,0);

    Scalar s2 = q2(0);
    Vector3 v2 = q2.block<3,1>(1,0);

    Vector3 v = v1.cross(v2) + s1 * v2 + s2 * v1;
    Scalar s = s1*s2 - v1.dot(v2);

    q << s, v;
    return q;
}

MEKF::Quaternion MEKF::quat_error_mult(const MEKF::Vector3 &a)
{
    Quaternion dq;

    // Factor
    const Scalar f = 1;

    // Scalar Part
    Scalar eta = 1 - a.squaredNorm() / 8;

    // Vector Part
    Vector3 nu = a / 2;

    // Error Parametrisation
    dq << eta, nu;
    dq *= f;

    return dq;
}

MEKF::MEKF(const SystemState &x0, const StateCov &P0) : P(P0)
{
    // initialize estimator state
    x.block<9,1>(0,0) = x0.block<9,1>(0,0);
    x.block<9,1>(9,0).setZero();

    // initialize reference quaternion
    qref << x0.block<4,1>(9,0);

    I.setIdentity();
}

void MEKF::predict(const Control &u)
{
    Eigen::Matrix<Scalar, nx, nx> F;
    SystemState xs;

    xs = get_system_state();

    // numerical integration
    xs = f.integrate(xs, u);

    std::cout << "xs = \n" << xs << std::endl;

    F = f.propagation_matrix(x, u, qref);

    std::cout << "F = \n" << F << std::endl;
    std::cout << "Q = \n" << f.Q << std::endl;

    x << xs.block<9,1>(0,0), F.block<9,18>(9,0)*x;

    // Covariance prediction
    P = F * P * F.transpose() + f.Q;
}

void MEKF::correct(const Control &u, const Measurement &z)
{
    Measurement y;                      // innovation
    Eigen::Matrix<Scalar, nz, nx> H;    // jacobian of h
    Eigen::Matrix<Scalar, nx, nx> IKH;  // temporary matrix
    Vector3 a;                          // attitude error parametrisation
    Quaternion dq;                      // attitude error quaternion

    H = h.jacobian(x, u, qref);
    std::cout << "H = \n" << H << std::endl;
    S = H * P * H.transpose() + h.R;
    std::cout << "S = \n" << S << std::endl;
    std::cout << "R = \n" << h.R << std::endl;

    // efficiently compute: K = P * H.transpose() * S.inverse();
    K = S.llt().solve(H * P).transpose();
    std::cout << "K = \n" << K << std::endl;

    y = z - h(x, u, qref);
    IKH = (I - K * H);
    std::cout << "y = \n" << y << std::endl;

    // Measurement update
    x = x + K * y;
    P = IKH * P * IKH.transpose() + K * h.R * K.transpose();

    // Attitude error transfer to reference quaternion qref
    a << x.block<3,1>(9,0);
    dq = quat_error_mult(a);
    qref = quatmul(dq, qref);

    // enforce unit norm constraint
    qref /= qref.norm();
    std::cout << "qref = \n" << qref << std::endl;

    // reset MEKF error angle
    x.block<3,1>(9,0).setZero();
}

void MEKF::update(const Control &u, const Measurement &z)
{
    predict(u);
    correct(u, z);
}

} // namespace bifoiler
