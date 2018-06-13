#include <iostream>
#include <Eigen/Dense>

#include "MEKF.h"
#include "model.h"

using namespace bifoiler;

int main()
{
    MEKF::SystemState x0;
    MEKF::StateCov P0;

    x0 << 5, 0, 0,    // v0 in BRF
          0, 0, 0,    // w0
          0, 0, -0.2, // r0 in IRF (NED)
          1, 0, 0, 0; // q0
    std::cout << "x0 = \n" << x0 << "\n";

    P0.setZero();
    P0.diagonal() << 0.3, 0.8, 0.5, // v
                     0.5, 0.5, 0.5, // w
                     0.1, 0.1, 0.3, // r
                     0.8, 0.8, 0.8, // a
                     0.01, 0.01, 0.01, // bg
                     0.01, 0.01, 0.01; // ba

    std::cout << "P0 = \n" << P0 << "\n";

    std::cout << "estimator test\n";
    MEKF estimator(x0, P0);

    MEKF::SystemState xs;
    MEKF::EstimatorState xe;
    MEKF::StateCov P;

    xs = estimator.get_system_state();
    std::cout << "x = \n" << xs << "\n";

    MEKF::Control u;
    MEKF::Measurement z;


    u << 0, 0, 0, 0.8f; // Flaps, Aileron, Rudder, Thrust
    z << 6, 0, 0, // v_B
         0, 0, 0, // w_B
         0.1, 0, 0, // r_I
         0, 0, -9.81; // acc_B

    std::cout << "foo\n";
    estimator.update(u, z);

    xs = estimator.get_system_state();
    std::cout << "xs = \n" << xs << "\n";

    xe = estimator.get_estimator_state();
    std::cout << "xe = \n" << xe << "\n";

    P = estimator.get_state_covariance();
    std::cout << "P = \n" << P << "\n";

}
