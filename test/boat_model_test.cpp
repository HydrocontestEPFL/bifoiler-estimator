
#include <casadi/casadi.hpp>
#include <boat_model.h>

using namespace casadi;
using namespace bifoiler;

namespace bifoiler {
    extern DM x0, u0;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " path/to/config.yaml" << std::endl;
        return -1;
    }

    std::string config_file(argv[1]);
    auto prop = BoatProperties::Load(config_file);

    x0 = DM::vertcat({5, 0, 0, // v
                         0.1, 0.1, 0.1, // W
                         0, 0, -0.2231, // r
                         1, 0, 0, 0 // q
    });
    u0 = DM::vertcat({0.1, 0.1, 0.1});

    BoatDynamics boat(prop);

    Function dynamics = boat.getNumericDynamics();
    std::cout << dynamics(SXVector{x0, u0}) << "\n";

    Function jacobian = boat.getNumericJacobian();
    std::cout << jacobian(SXVector{x0, u0}) << "\n";

    Dict opt = {
        {"with_header", true},
    };
    dynamics.generate("codegen_dynamics", opt);
}
