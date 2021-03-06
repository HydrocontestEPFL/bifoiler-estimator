#include <iostream>
#include <fstream>
#include <string>

#include <casadi/casadi.hpp>
#include <boat_model.h>
#include <estimator_model.h>

using namespace casadi;
using namespace bifoiler;

#define DEG2RAD(x) ((x) * M_PI / 180)

void csv_write_line(std::ofstream &out, DM x)
{
    for (int i = 0; i < x.size1(); i++) {
        if (i != 0) {
            out << ",";
        }
        out << x(i);
    }
    out << "\n";
}

DM cvodes_solve(Function &cvodes_integrator, const DM &x0, const DM &u)
{
    DMDict out;
    try
    {
        DMDict args = {{"x0", x0}, {"p", u}};
        out = cvodes_integrator(args);
    }
    catch(std::exception &e)
    {
        std::cout << "At state x0 : " << x0 << "\n";
        std::cout << "At control u:" << u << "\n";
        std::cout << "CVODES exception " << e.what() << "\n";
    }

    return out["xf"];
}

DM rk4_solve(Function &integrator, const DM &x0, const DM &u, const double h)
{
    DMVector arg = DMVector{x0,u,h};
    auto out = integrator(arg);
    return out[0];
}

DM h_func(Function h, DM xs, DM u)
{
    DM qr = xs(Slice(9,13));
    DM xe = DM::vertcat({
        xs(Slice(0,9)),
        DM::zeros(3,1), // a attitude error parametrization
        DM::zeros(3,1), // bg gyro bias
        DM::zeros(3,1)  // ba accelerometer bias
    });

    DMVector z = h(DMVector{xe, u, qr});
    return z[0];
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " config.yaml [RK4|CVODES]" << std::endl;
        return -1;
    }

    std::string integrator_name = "RK4";
    if (argc >= 3) {
        integrator_name = argv[2];
    }

    std::ofstream sim_x("sim_x.csv");
    std::ofstream sim_u("sim_u.csv");
    std::ofstream sim_z("sim_z.csv");

    sim_x << "vx,vy,vz,wx,wy,wz,rx,ry,rz,q0,q1,q2,q3\n";
    sim_u << "flaps,aileron,rudder,thrust\n";
    sim_z << "vx,vy,vz,wx,wy,wz,rx,ry,rz,ax,ay,az\n";


    std::cout << "load config file: " << argv[1] << std::endl;
    std::string config_file(argv[1]);
    BoatProperties prop = BoatProperties::Load(config_file);

    BoatDynamics boat_model(prop);
    EstimatorModel estimator(boat_model, prop);

    SX dynamics = boat_model.getSymbolicDynamics();
    SX state = boat_model.getSymbolicState();
    SX control = boat_model.getSymbolicControl();   // Control: Flaps, Ailerons, Rudder

    Function output_map = estimator.h_func;

    const double h = prop.estimator.t_samp;
    SXDict ode = {{"x", state}, {"p", control}, {"ode", dynamics}};
    Dict opts = {{"tf", h}};
    Function CVODES_INT = integrator("CVODES_INT", "cvodes", ode, opts);
    Function RK4_INT = boat_model.getNumericIntegrator();

    const double SIM_TIME = 10; // [s]
    DM x = DM::vertcat({
        5, 0, 0,    // v0 [m/s] in BRF
        0, 0, 0,    // w0 [rad/s]
        0, 0, 0, // r0 [m] in IRF (NED)
        1, 0, 0, 0  // q0
    });

    std::cout << "start simulation" << std::endl;
    for (double t = 0; t < SIM_TIME; t += h) {
        DM u = DM::vertcat({
            0, // Flaps
            0, // Aileron
            0, // Rudder
            // DEG2RAD(3)*sin(2*M_PI*0.5*t), // Flaps
            // DEG2RAD(3)*sin(2*M_PI*0.5*t), // Aileron
            // DEG2RAD(3)*sin(2*M_PI*0.1*t), // Rudder
            0.8 // Thrust
        });

        DM z = h_func(output_map, x, u);

        // std::cout << x << std::endl;
        // std::cout << z << std::endl;

        csv_write_line(sim_x, x);
        csv_write_line(sim_u, u);
        csv_write_line(sim_z, z);

        if (integrator_name == "CVODES") {
            x = cvodes_solve(CVODES_INT, x, u);
        } else if (integrator_name == "RK4") { // RK4
            x = rk4_solve(RK4_INT, x, u, h);
        } else {
            std::cerr << "unknown integrator: " << integrator_name << std::endl;
            return -1;
        }
    }

    std::cout << "DONE" << std::endl;
}
