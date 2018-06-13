#include <casadi/casadi.hpp>

#include "estimator_model.h"
#include "boat_math.h"

using namespace casadi;
using namespace boat_math;

namespace bifoiler {

void quat_error_mult(const SX &a, const SX &qr, SX &q, SX &eta)
{
    // Attitude error transfer to reference quaternion q
    const double f = 1;
    eta = 1 - SX::dot(a, a) / 8;
    SX nu = a / 2;
    SX dq = f * SX::vertcat({eta, nu});

    q = quatmul(dq, qr);
    q /= SX::norm_2(q); // enforce unit norm constraint
}

EstimatorModel::EstimatorModel(BoatDynamics &dynamics, const BoatProperties &prop)
{
    const double t_samp = prop.estimator.t_samp;

    SX a = SX::sym("a", 3);      // Attitude error parametrisation
    SX bg = SX::sym("bg", 3);    // IMU accelerometer bias
    SX ba = SX::sym("ba", 3);    // IMU gyroscope bias
    SX qr = SX::sym("qr", 4);    // reference, true of last iteration

    SX xs = dynamics.getSymbolicState();
    SX u = dynamics.getSymbolicControl();   // Control: Flaps, Ailerons, Rudder

    // Erroneous Quaternion
    SX q, eta;
    quat_error_mult(a, qr, q, eta);

    SX v = xs(Slice(0,3));
    SX W = xs(Slice(3,6));
    SX r = xs(Slice(6,9));
    SX q_xs = xs(Slice(9,13));
    SX xe = SX::vertcat({v, W, r, a, bg, ba}); // estimator state

    // (M)EKF dynamics/ Composition of total Jacobian

    // Model Jacobian to get linearized model
    SX A_sys = substitute(dynamics.getSymbolicJacobian(), q_xs, q);

    // df/da = df/dq*dq/da; Jq := dq/da
    SX Jq = SX::jacobian(q, a);

    // Matrix valued chain rule
    SX Aa = SX::mtimes(A_sys(Slice(0, 9), Slice(9, 13)), Jq);

    // Attitude Dynamics (from Nokland)
    SX f_att_nl = SX::vertcat({
        0.5 * SX::mtimes(SX::eye(3) * eta + skew_mat(a), W - bg),
        SX::zeros(3,1)
    });

    // Attitude Dynamics (from Markley)
    SX f_att_mk = SX::vertcat({
        -bg - SX::mtimes(skew_mat(W), a),
        SX::zeros(3,1)
    });

    SX Jatt = SX::jacobian(f_att_mk, xe);

    SX A = SX::vertcat({
        SX::horzcat({A_sys(Slice(0, 9), Slice(0, 9)), Aa, SX::zeros(9, 6)}),
        Jatt,
        SX::zeros(3, xe.size1())
    });
    Function A_func = Function("A_func", {xe, u, qr}, {F});

    // Discretization
    SX F = SX::eye(xe.size1()) + t_samp*A; // Euler
    Function F_func = Function("F_func", {xe, u, qr}, {F});

    this->F = A;
    this->A_func = A_func;
    this->F = F;
    this->F_func = F_func;

    // Outputmap and Sensitivity
    SX g = SX::vertcat({0, 0, -prop.env.g}); // in BRF: z = up
    SX r_ant = SX::vertcat({
        prop.sensor.r_ant[0],
        prop.sensor.r_ant[1],
        prop.sensor.r_ant[2]
    });
    SX g_BRF = quatrot(q, g);

    // Velocity is in BRF
    SX y_vgns = quatrot_inverse(q, v) + quatrot_inverse(q, SX::mtimes(skew_mat(W), r_ant));
    SX y_wgyro = W + bg;
    SX y_rgns = r + quatrot_inverse(q, r_ant);

    //y_acc = v_dot_BRF - quatrot(q,[0; g]) + ba;
    SX y_acc = -quatrot(q, SX::vertcat({0,g})) + ba;

    // Complete Output-map h(x) and Output Sensitivity H
    SX h = SX::vertcat({y_vgns, y_wgyro, y_rgns, y_acc});
    Function h_func = Function("h_func", {xe, u, qr}, {h});
    SX H = SX::jacobian(h, xe);
    Function H_func = Function("H_func", {xe, u, qr}, {H});

    this->h = h;
    this->h_func = h_func;
    this->H = H;
    this->H_func = H_func;
}

} // namespace bifoiler
