#ifndef MODEL_H
#define MODEL_H

#include <Eigen/Dense>

namespace bifoiler {

struct Dynamics {
    using Scalar = double;
    using SystemState = Eigen::Matrix<Scalar, 13, 1>;
    using EstimatorState = Eigen::Matrix<Scalar, 18, 1>;
    using Control = Eigen::Matrix<Scalar, 4, 1>;
    using Quaternion = Eigen::Matrix<Scalar, 4, 1>; // [w, x, y, z];
    enum {
        nx = EstimatorState::RowsAtCompileTime,
        nxs = SystemState::RowsAtCompileTime,
        nu = Control::RowsAtCompileTime,
    };

    Eigen::Matrix<Scalar, nx, nx> Q;

    Dynamics();
    SystemState integrate(const SystemState &x, const Control &u);
    Eigen::Matrix<Scalar, nx, nx> propagation_matrix(const EstimatorState &x,
                                                     const Control &u,
                                                     const Quaternion &qref);
};

struct Observation {
    using Scalar = double;
    using EstimatorState = Eigen::Matrix<Scalar, 18, 1>;
    using Control = Eigen::Matrix<Scalar, 4, 1>;
    using Measurement = Eigen::Matrix<Scalar, 12, 1>;
    using Quaternion = Eigen::Matrix<Scalar, 4, 1>; // [w, x, y, z];
    enum {
        nx = EstimatorState::RowsAtCompileTime,
        nu = Control::RowsAtCompileTime,
        nz = Measurement::RowsAtCompileTime,
    };

    Eigen::Matrix<Scalar, nz, nz> R;

    Observation();
    Measurement operator()(const EstimatorState &x,
                           const Control &u,
                           const Quaternion &qref);
    Eigen::Matrix<Scalar, nz, nx> jacobian(const EstimatorState &x,
                                           const Control &u,
                                           const Quaternion &qref);
};

} // namespace bifoiler

#endif /* MODEL_H */
