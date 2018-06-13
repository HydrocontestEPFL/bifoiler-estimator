#include <unistd.h>

#include <casadi/casadi.hpp>
#include <boat_model.h>
#include <estimator_model.h>

using namespace casadi;
using namespace bifoiler;

static void generate(Function func, const char *dir, const char *name, Dict &opt)
{
    std::cout << "generate: " << dir << "/" << name << ".c\n";
    func.generate(name, opt);
}

int main(int argc, char *argv[])
{
    char *dir = NULL;
    const char *outdir = ".";

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


        outdir = argv[2];
        dir = getcwd(buf, sizeof(buf));
        err = chdir(outdir);

        if (err != 0) {
            std::cout << "error can't chdir to directory: " << argv[2] << std::endl;
            std::cout << "does it exist?" << std::endl;
            return -1;
        }
    }

    BoatDynamics boat_model(prop);
    EstimatorModel estimator(boat_model, prop);

    Dict opt = {
        {"with_header", true},
    };

    generate(boat_model.getNumericIntegrator(), outdir, "integrator", opt);
    generate(estimator.A_func, outdir, "estimator_jacobian", opt);
    generate(estimator.F_func, outdir, "propagation_matrix", opt);
    generate(estimator.h_func, outdir, "output_map", opt);
    generate(estimator.H_func, outdir, "output_map_jacobian", opt);

    // return to old working directory
    if (dir != NULL) {
        chdir(dir);
    }

    return 0;
}
