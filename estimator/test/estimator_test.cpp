#include <iostream>
#include <Eigen/Dense>

#include "MEKF.h"
#include "model.h"

using namespace bifoiler;

using KF = MEKF<Dynamics, Observation>;

int main()
{
    KF::SystemState x0;
    KF::StateCov P0;

    x0 << 5, 0, 0,    // v0 in BRF
          0, 0, 0,    // w0
          0, 0, -0.2, // r0 in IRF (NED)
          1, 0, 0, 0; // q0

    P0.setZero();
    P0.diagonal() << 0.3, 0.8, 0.5, // v
                     0.5, 0.5, 0.5, // w
                     0.1, 0.1, 0.3, // r
                     0.8, 0.8, 0.8, // a
                     0.01, 0.01, 0.01, // bg
                     0.01, 0.01, 0.01; // ba

    std::cout << "estimator test\n";
    KF estimator(x0, P0);

    KF::SystemState xs;
    KF::EstimatorState xe;
    KF::StateCov P;

    xs = estimator.get_system_state();
    std::cout << "x = " << xs << "\n";

    KF::Control u;
    KF::Measurement z;

    u << 0, 0, 0, 0.8f; // Flaps, Aileron, Rudder, Thrust
    z << 0, 0, 0, 0.8f; //

    estimator.update(u, z);

    xs = estimator.get_system_state();
    std::cout << "xs = " << xs << "\n";

    xe = estimator.get_estimator_state();
    std::cout << "xe = " << xe << "\n";

    P = estimator.get_state_covariance();
    std::cout << "P = " << P << "\n";

}
