#include <unistd.h>

#include <casadi/casadi.hpp>
#include <boat_model.h>
#include <estimator_model.h>

using namespace casadi;
using namespace bifoiler;

int main(int argc, char *argv[])
{
    char *dir = NULL;

    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " config.yaml [codegen/dir]" << std::endl;
        return -1;
    }

    std::string config_file(argv[1]);
    std::cout << "load config file: " << config_file << std::endl;
    BoatProperties prop = BoatProperties::Load(config_file);

    // change working directory
    if (argc >= 3) {
        char buf[1000];
        int err;

        dir = getcwd(buf, sizeof(buf));
        err = chdir(argv[2]);

        if (err != 0) {
            std::cout << "error with codegen directory: " << argv[2] << std::endl;
            std::cout << "does it exist?" << std::endl;
            return -1;
        }
    }

    BoatDynamics boat_model(prop);
    EstimatorModel estimator(boat_model, prop);

    Dict opt = {
        {"with_header", true},
        {"verbose", true},
    };
    std::cout << "generate: integrator.c\n";
    boat_model.getNumericIntegrator().generate("integrator", opt);

    std::cout << "F =\n" << estimator.F << "\n";
    std::cout << "generate: F_func.c\n";
    estimator.F_func.generate("F_func", opt);

    // std::cout << "h =\n" << estimator.h << "\n";
    std::cout << "generate: h_func.c\n";
    estimator.h_func.generate("h_func", opt);

    // std::cout << "H =\n" << estimator.H << "\n";
    std::cout << "generate: H_func.c\n";
    estimator.H_func.generate("H_func", opt);


    // return to old working directory
    if (dir != NULL) {
        chdir(dir);
    }

    return 0;
}
