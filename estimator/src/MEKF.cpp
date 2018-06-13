#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "MEKF.h"

namespace bifoiler {

template <typename Dynamics, typename Observation>
typename MEKF<Dynamics, Observation>::SystemState MEKF<Dynamics, Observation>::get_system_state()
{
    SystemState xs;
    xs << x.block<9,1>(0,0), qref;
    return xs;
}

template <typename Dynamics, typename Observation>
typename MEKF<Dynamics, Observation>::Quaternion MEKF<Dynamics, Observation>::quatmul(const MEKF<Dynamics, Observation>::Quaternion &q1, const MEKF<Dynamics, Observation>::Quaternion &q2)
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

template <typename Dynamics, typename Observation>
typename MEKF<Dynamics, Observation>::Quaternion MEKF<Dynamics, Observation>::quat_error_mult(const MEKF<Dynamics, Observation>::Vector3 &a)
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

template <typename Dynamics, typename Observation>
MEKF<Dynamics, Observation>::MEKF(const SystemState &x0, const StateCov &P0) : P(P0)
{
    x << x0.block<9,1>(0,0);
    x.block<9,1>(9,0).setZero();
    I.setIdentity();
}

template <typename Dynamics, typename Observation>
void MEKF<Dynamics, Observation>::predict(const Control &u)
{
    Eigen::Matrix<Scalar, nx, nx> F;
    SystemState xs;

    xs = get_system_state();

    // numerical integration
    xs = f.integrate(xs, u);

    F = f.propagation_matrix(x, u, qref);

    x << xs.block<9,1>(0,0), F.block<9,18>(9,0)*x;

    // Covariance prediction
    P = F * P * F.transpose() + f.Q;
}

template <typename Dynamics, typename Observation>
void MEKF<Dynamics, Observation>::correct(const Control &u, const Measurement &z)
{
    Measurement y;                      // innovation
    Eigen::Matrix<Scalar, nz, nx> H;    // jacobian of h
    Eigen::Matrix<Scalar, nx, nx> IKH;  // temporary matrix
    Vector3 a;                          // attitude error parametrisation
    Quaternion dq;                      // attitude error quaternion

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
    a << x.block<3,1>(9,0);
    dq = quat_error_mult(a);
    qref = quatmul(dq, qref);

    // enforce unit norm constraint
    qref /= qref.norm();

    // reset MEKF error angle
    x.block<3,1>(9,0).setZero();
}

template <typename Dynamics, typename Observation>
void MEKF<Dynamics, Observation>::update(const Control &u, const Measurement &z)
{
    predict(u);
    correct(u, z);
}

} // namespace bifoiler
