#include <iostream>
#include "codegen_dynamics.h"
#include <Eigen/Dense>

// dynamics {13, 3} -> {13}

struct Dynamics {
    using State = Eigen::Matrix<double, 13, 1>;
    using Control = Eigen::Matrix<double, 3, 1>;
    enum {
        nx = State::RowsAtCompileTime,
        nu = Control::RowsAtCompileTime,
    };

    State operator()(const State &x, const Control &u)
    {
        double _x[nx], _u[nu], _xd[nx];

        State::Map(_x) = x;
        Control::Map(_u) = u;

        // Setup function arguments
        // TODO: verify correct usage
        const double *arg[2] = {_x, _u};
        double *res[1] = {_xd};
        int* iw = NULL;
        double* w = NULL;
        void* mem = NULL;

        dynamics(&arg[0], &res[0], iw, w, mem);

        return State::Map(_xd);
    }
};

using State = Dynamics::State;
using Control = Dynamics::Control;

int main()
{
    State x, xd;
    Control u;
    Dynamics f;

    x << 5, 0, 0, // v
         0.1, 0.1, 0.1, // W
         0, 0, -0.2231, // r
         1, 0, 0, 0; // q

    u << 0.1, 0.1, 0.1;

    xd = f(x, u);

    std::cout << xd.transpose() << "\n";
}
