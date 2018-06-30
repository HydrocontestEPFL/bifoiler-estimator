
#include <codegen/integrator.h>
#include <codegen/propagation_matrix.h>
#include <codegen/output_map.h>
#include <codegen/output_map_jacobian.h>

#include "model.h"
#include "covariance.h"

namespace bifoiler {

using Scalar = Dynamics::Scalar;

enum {
    nx = Dynamics::nx,
    nxs = Dynamics::nxs,
    nu = Dynamics::nu,
    nz = Observation::nz,
};

Dynamics::Dynamics()
{
    Q = process_noise_covariance<Scalar>();
}

Dynamics::SystemState Dynamics::integrate(const SystemState &x, const Control &u)
{
    const unsigned N = 10;
    Scalar _x0[nxs], _u[nu], _dT[1], _x1[nxs];

    SystemState::Map(_x0) = x;
    Control::Map(_u) = u;
    _dT[0] = 0.02 / N; // TODO: parametrize

    // Setup function arguments
    const Scalar *arg[3] = {_x0, _u, _dT};
    Scalar *res[1] = {_x1};

    for (unsigned i = 0; i < N; i++) {
        // call to CasADi generated C function
        // integrator {13, 3, 1} -> {13}
        integrator(&arg[0], &res[0], NULL, NULL, NULL);

        // copy output to input
        SystemState::Map(_x0) = SystemState::Map(_x1);
    }

    return SystemState::Map(_x1);
}

Eigen::Matrix<Scalar, nx, nx> Dynamics::propagation_matrix(const EstimatorState &x,
                                                 const Control &u,
                                                 const Quaternion &qref)
{
    Scalar _x[nx], _u[nu], _qref[4], _F[nx*nx];

    EstimatorState::Map(_x) = x;
    Control::Map(_u) = u;
    Quaternion::Map(_qref) = qref;

    // Setup function arguments
    const Scalar *arg[3] = {_x, _u, _qref};
    Scalar *res[1] = {_F};

    // call to CasADi generated C function
    // F_func {18, 3, 4} -> {18}
    F_func(&arg[0], &res[0], NULL, NULL, NULL);

    return Eigen::Matrix<Scalar, nx, nx>::Map(_F);
}

Observation::Observation()
{
    R = measurement_noise_covariance<Scalar>();
}

Observation::Measurement Observation::operator()(const EstimatorState &x,
                       const Control &u,
                       const Quaternion &qref)
{
    Scalar _x[nx], _u[nu], _qref[4], _z[nz];

    EstimatorState::Map(_x) = x;
    Control::Map(_u) = u;
    Quaternion::Map(_qref) = qref;

    // Setup function arguments
    const Scalar *arg[3] = {_x, _u, _qref};
    Scalar *res[1] = {_z};

    // call to CasADi generated C function
    // h_func {18, 3, 4} -> {12}
    h_func(&arg[0], &res[0], NULL, NULL, NULL);

    return Measurement::Map(_z);
}

Eigen::Matrix<Scalar, nz, nx> Observation::jacobian(const EstimatorState &x,
                                       const Control &u,
                                       const Quaternion &qref)
{
    Scalar _x[nx], _u[nu], _qref[4], _H[nz*nx];

    EstimatorState::Map(_x) = x;
    Control::Map(_u) = u;
    Quaternion::Map(_qref) = qref;

    // Setup function arguments
    const Scalar *arg[3] = {_x, _u, _qref};
    Scalar *res[1] = {_H};

    // call to CasADi generated C function
    // h_func {18, 3, 4} -> {12x18}
    H_func(&arg[0], &res[0], NULL, NULL, NULL);

    return Eigen::Matrix<Scalar, nz, nx>::Map(_H);
}

} // namespace bifoiler
