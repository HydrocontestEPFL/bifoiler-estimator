#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <Eigen/Dense>

#include "MEKF.h"
#include "model.h"
#include "covariance.h"

using namespace bifoiler;
using namespace Eigen;

MEKF::Measurement add_noise(MEKF::Measurement z)
{
    static std::default_random_engine generator;
    static std::normal_distribution<double> noise(0, 1);

    auto Q = measurement_noise_covariance<double>();
    MEKF::Measurement variance = Q.diagonal();

    for (int i = 0; i < z.size(); i++) {
        z(i) += sqrt(variance(i)) * noise(generator); // additive gaussian noise
    }
    return z;
}

template <typename vector>
int csv_parse_line(vector &out, std::ifstream &in)
{
    std::string line;
    if (!std::getline(in, line)) {
        return -1;
    }

    std::stringstream ss(line);

    unsigned i = 0;
    double d;
    while (ss >> d && i < out.size())
    {
        out(i) = d;

        if (ss.peek() == ',') {
            ss.ignore();
        }
        i++;
    }
    if (i != out.size()) {
        std::cout << "CSV parse error: " << line << std::endl;
        return -1;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " path/to/csv/data" << std::endl;
        return -1;
    }

    std::string sim_csv_path(argv[1]);
    std::ifstream sim_x(sim_csv_path + "/sim_x.csv");
    std::ifstream sim_u(sim_csv_path + "/sim_u.csv");
    std::ifstream sim_z(sim_csv_path + "/sim_z.csv");

    // quick hack to discard CSV data header
    std::string line;
    std::getline(sim_x, line);
    std::getline(sim_u, line);
    std::getline(sim_z, line);

    MEKF::SystemState xs, xs_sim;
    MEKF::EstimatorState xe;
    MEKF::StateCov P;
    MEKF::Control u;
    MEKF::Measurement z;

    MEKF::SystemState x0;
    MEKF::StateCov P0;

    x0 << 5, 0, 0,    // v0 in BRF
          0, 0, 0,    // w0
          0, 0, 0,    // r0 in IRF (NED)
          1, 0, 0, 0; // q0

    P0.setZero();
    P0.diagonal() << 0.3, 0.8, 0.5, // v
                     0.5, 0.5, 0.5, // w
                     0.1, 0.1, 0.3, // r
                     0.8, 0.8, 0.8, // a
                     0.01, 0.01, 0.01, // bg
                     0.01, 0.01, 0.01; // ba


    bool run = true;
    run &= (csv_parse_line<MEKF::SystemState>(x0, sim_x) == 0);

    std::cout << "init MEKF" << std::endl;
    MEKF estimator(x0, P0);

    while (run) {
        run &= (csv_parse_line<MEKF::Control>(u, sim_u) == 0);
        run &= (csv_parse_line<MEKF::Measurement>(z, sim_z) == 0);
        // std::cout << "u = \n" << u.transpose() << "\n";
        // std::cout << "z = \n" << z.transpose() << "\n";

        z = add_noise(z);

        estimator.update(u, z);

        xs = estimator.get_system_state();
        std::cout << "xs  = " << xs.transpose() << "\n";

        run &= (csv_parse_line<MEKF::SystemState>(xs_sim, sim_x) == 0);
        std::cout << "sim = " << xs.transpose() << "\n";

        for (int i = 0; i < xs.size(); i++) {
            if (isnan(xs(i))) {
                std::cout << "Nan" << std::endl;
                run = false;
                break;
            }
        }
    }
    std::cout << "DONE" << std::endl;
}
