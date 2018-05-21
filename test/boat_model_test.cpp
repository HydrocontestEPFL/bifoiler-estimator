
#include <casadi/casadi.hpp>
#include <boat_model.h>

using namespace casadi;
using namespace bifoiler;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " path/to/config.yaml" << std::endl;
        return -1;
    }

    std::string config_file(argv[1]);
    auto prop = BoatProperties::Load(config_file);

    SX v    = SX::sym("v", 3);   // [m/s] translational velocity of the CoG in BRF
    SX W    = SX::sym("W", 3);   // [rad/s] angular velocities of boat in BRF
    SX r    = SX::sym("r", 3);   // [m] position of the CoG in IRF
    SX q_BI = SX::sym("q", 4);   // unit quaternion for transformation from IRF to BRF
    SX state    = SX::vertcat({v, W, r, q_BI});

    SX dF   = SX::sym("dF");     // Flaps
    SX dA   = SX::sym("dA");     // Aileron deflection [reserved, but not used] [rad]
    SX dR   = SX::sym("dR");     // Rudder deflection [rad]
    SX dE   = SX::sym("dE");     // Elevator deflection [positive down] [rad]
    SX control  = SX::vertcat({dF, dA, dR, dE});

    SX Fhbrf, Mhbrf, aoa, ssa;

    BoatDynamics::Hydrodynamics(state, control, prop, Fhbrf, Mhbrf, aoa, ssa);

    std::cout << Fhbrf << "\n" << Mhbrf << "\n" << aoa << "\n" << ssa << "\n";
}
