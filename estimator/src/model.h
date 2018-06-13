#ifndef MODEL_H
#define MODEL_H

#include <codegen/integrator.h>
#include <codegen/propagation_matrix.h>
#include <codegen/output_map.h>
#include <codegen/output_map_jacobian.c.h>
#include <Eigen/Dense>
#include "covariance.h"

namespace bifoiler {

using Quaternion = Eigen::Matrix<Scalar, 4, 1>; // [w, x, y, z];

struct Dynamics {
    using Scalar = double;
    using SystemState = Eigen::Matrix<Scalar, 13, 1>;
    using EstimatorState = Eigen::Matrix<Scalar, 13, 1>;
    using Control = Eigen::Matrix<Scalar, 3, 1>;
    enum {
        nx = EstimatorState::RowsAtCompileTime,
        nxs = SystemState::RowsAtCompileTime,
        nu = Control::RowsAtCompileTime,
    };

    Eigen::Matrix<Scalar, nx, nx> Q;

    Dynamics()
    {
        Q = process_noise_covariance<Scalar>();
    }

    SystemState integrate(const SystemState &x, const Control &u)
    {
        Scalar _x0[nxs], _u[nu], _x1[nxs], _dT[1];

        SystemState::Map(_x0) = x;
        Control::Map(_u) = u;
        _dT[0] = 0.02; // TODO: parametrize

        // Setup function arguments
        const Scalar *arg[2] = {_x0, _u, _dT};
        Scalar *res[1] = {_x1};

        // call to CasADi generated C function
        // integrator {13, 3, 1} -> {13}
        integrator(&arg[0], &res[0], NULL, NULL, NULL);

        return SystemState::Map(_x1);
    }

    Eigen::Matrix<Scalar, nx, nx> propagation_matrix(const EstimatorState &x,
                                                     const Control &u,
                                                     const Quaternion &qref)
    {
        Scalar _x[nx], _u[nu], _qref[4], _F[nx*nx];

        SystemState::Map(_x0) = x;
        Control::Map(_u) = u;
        Quaternion::Map(_qref) = qref;

        // Setup function arguments
        Scalar *arg[2] = {_x0, _u, _qref};
        Scalar *res[1] = {_F};

        // call to CasADi generated C function
        // F_func {18, 3, 4} -> {18}
        F_func(&arg[0], &res[0], NULL, NULL, NULL);

        return Eigen::Matrix<Scalar, nx, nx>::Map(_F);
    }
};

struct Observation {
    using Scalar = double;
    using EstimatorState = Eigen::Matrix<Scalar, 13, 1>;
    using Control = Eigen::Matrix<Scalar, 3, 1>;
    using Measurement = Eigen::Matrix<Scalar, 12, 1>;
    enum {
        nx = EstimatorState::RowsAtCompileTime,
        nu = Control::RowsAtCompileTime,
        nz = Measurement::RowsAtCompileTime,
    };

    Eigen::Matrix<Scalar, nz, nz> R;

    Observation()
    {
        R = measurement_noise_covariance<Scalar>();
    }

    Measurement operator()(const EstimatorState &x,
                           const Control &u,
                           const Quaternion &qref)
    {
        Scalar _x[nx], _u[nu], _qref[4], _z[nz];

        EstimatorState::Map(_x) = x;
        Control::Map(_u) = u;
        Quaternion::Map(_qref) = qref;

        // Setup function arguments
        const Scalar *arg[2] = {_x, _u, _qref};
        Scalar *res[1] = {_z};

        // call to CasADi generated C function
        // h_func {18, 3, 4} -> {12}
        h_func(&arg[0], &res[0], NULL, NULL, NULL);

        return Measurement::Map(_z);
    }

    Eigen::Matrix<Scalar, nz, nx> jacobian(const EstimatorState &x,
                                           const Control &u,
                                           const Quaternion &qref)
    {
        Scalar _x[nx], _u[nu], _qref[4], _H[nz*nx];

        EstimatorState::Map(_x0) = x;
        Control::Map(_u) = u;
        Control::Map(_qref) = qref;

        // Setup function arguments
        Scalar *arg[2] = {_x0, _u, _qref};
        Scalar *res[1] = {_H};

        // call to CasADi generated C function
        // h_func {18, 3, 4} -> {12x18}
        H_func(&arg[0], &res[0], NULL, NULL, NULL);

        return Eigen::Matrix<Scalar, nz, nx>::Map(_H);
    }
};

} // namespace bifoiler

#endif /* MODEL_H */
