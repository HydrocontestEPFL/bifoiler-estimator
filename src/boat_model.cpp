#include <iostream>
#include <math.h>

#include <casadi/casadi.hpp>

#include "boat_model.h"
#include "boat_properties.h"

using namespace casadi;

namespace bifoiler {

SX quatmul(const SX &q1, const SX &q2)
{
    SX s1 = q1(0);
    SX v1 = q1(Slice(1,4));

    SX s2 = q2(0);
    SX v2 = q2(Slice(1,4));

    SX q = SX::vertcat({s1*s2 - SX::dot(v1,v2), SX::cross(v1,v2) + s1 * v2 + s2 * v1});
    return q;
}

SX quatconj(const SX &q)
{
    return SX::vertcat({q(0), -q(1), -q(2), -q(3)});
}

SX quatinv(const SX &q)
{
    // TODO: should be divided by squared norm
    return quatconj(q) / SX::norm_2(q);

    // SX norm_squared = SX::dot(q,q);
    // return quatconj(q) / norm_squared;
}

// TODO: verify correctness (multiplication order)
SX quatrot(const SX &q, const SX &r)
{
    // SX qinv = quatinv(q);
    SX qinv = quatconj(q);  // q_inv = q_conj for unit quaternions
    SX qr   = SX::vertcat({0, r});
    SX qrr  = quatmul(quatmul(qinv, qr),q);
    return qrr(Slice(1,4));
}

SX quatrot_inverse(const SX &q, const SX &r)
{
    // SX qinv = quatinv(q);
    SX qinv = quatconj(q);  // q_inv = q_conj for unit quaternions
    SX qr   = SX::vertcat({0, r});
    SX qrr  = quatmul(quatmul(q, qr),qinv);
    return qrr(Slice(1,4));
}

SX rk4_symbolic(const SX &x,
                const SX &u,
                Function &func,
                const SX &h)
{
    SXVector res = func(SXVector{x, u});
    SX k1 = res[0];
    res = func(SXVector{x + 0.5 * h * k1, u});
    SX k2 = res[0];
    res = func(SXVector{x + 0.5 * h * k2, u});
    SX k3 = res[0];
    res = func(SXVector{x + h * k3, u});
    SX k4 = res[0];

    return x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
}

SX skew_mat(const SX &v)
{
    return SX::vertcat({
        SX::horzcat({0, -v(2), v(1)}),
        SX::horzcat({v(2), 0, -v(0)}),
        SX::horzcat({-v(1), v(0), 0})
    });
}

// TODO: debug only, remove later
DM x0, u0;
static void print_expr(const char *name, SX &sym, SX &x, SX &u)
{
    std::cout << name << "\n";
    std::cout << sym << "\n";
    Function f = Function(name,{x, u},{sym});
    std::cout << f(SXVector{x0, u0}) << "\n";
}

// TODO: debug only, remove later
static void print_jacobian(const char *name, SX &sym, SX &x, SX &u)
{
    SX jac = SX::jacobian(sym, x);
    print_expr(name, jac, x, u);
}

void BoatDynamics::Hydrodynamics(const SX &state, const SX &control, const BoatProperties &prop, SX &Fhbrf, SX &Mhbrf, SX &aoa, SX &ssa)
{
    double rho = prop.env.rho_sh2o;
    double AR  = prop.foils.ARff;        // [-], Aspect Ratio (Flügelstreckung)
    double S   = prop.foils.areaff;      // [m^2], Frontfoil Wing area
    double c   = prop.foils.mac;         // [m], Mean aerodynamic chord
    double b   = prop.foils.wingspanff;

    // Static Hydrodynamic Coefficient
    // All characteristics assumed linear
    double CL0         = prop.hydrodynamic.CL0;
    double CLa_tot     = prop.hydrodynamic.CLa_total;
    double e_o         = prop.hydrodynamic.e_oswald;
    double dw          = CLa_tot / (M_PI * e_o * AR); // downwash acting at the tail []
    double CD0_tot     = prop.hydrodynamic.CD0_total;
    double CYb         = prop.hydrodynamic.CYb;
    double Cm0         = prop.hydrodynamic.Cm0;
    double Cma         = prop.hydrodynamic.Cma;
    double Cn0         = prop.hydrodynamic.Cn0;
    double Cnb         = prop.hydrodynamic.Cnb;
    double Cl0         = prop.hydrodynamic.Cl0;
    double Clb         = prop.hydrodynamic.Clb;
    double CLq         = prop.hydrodynamic.CLq;
    double Cmq         = prop.hydrodynamic.Cmq;
    double CYr         = prop.hydrodynamic.CYr;
    double Cnr         = prop.hydrodynamic.Cnr;
    double Clr         = prop.hydrodynamic.Clr;
    double CYp         = prop.hydrodynamic.CYp;
    double Clp         = prop.hydrodynamic.Clp;
    double Cnp         = prop.hydrodynamic.Cnp;

    // Control Derivatives
    double CXdf        = prop.hydrodynamic.CXdf;
    double CYdr        = prop.hydrodynamic.CYdr;
    double CZde        = 0; // prop.hydrodynamic.CZde;
    double CZdf        = prop.hydrodynamic.CZdf;
    double Clda        = prop.hydrodynamic.CLda;
    double Cldr        = prop.hydrodynamic.CLdr;
    double Cmde        = 0; // prop.hydrodynamic.CMde;
    double Cmdf        = prop.hydrodynamic.CMdf;
    double Cnda        = prop.hydrodynamic.CNda;
    double Cndr        = prop.hydrodynamic.CNdr;

    // State variables definition
    SX v = state(Slice(0,3)); // Translational Velocity/Attitude towards Velocity Reference Frame
    SX W = state(Slice(3,6)); // Rotational Rate/Attitude Derivative

    // Control variables definition
    SX dF = control(0); // Flaps
    SX dA = control(1); // Aileron
    SX dR = control(2); // Rudder
    // SX dE = control(3); // Elevator

    // Variables and Orientation
    const double v1_eps   = 0.0001;    // Enable v(1)_0 = 0 w/o getting NaN
    SX V = SX::norm_2(v)+ v1_eps;        // Absolute Value of Velocity

    ssa = asin(v(1) / V);               // side slip angle [rad]
    aoa = atan2(v(2) , v(0) + v1_eps);  // angle of attack definition [rad]

    SX dyn_press = 0.5 * rho * pow(V,2);  // dynamic pressure

    SX q_aoa = SX::vertcat({cos(aoa/2), 0, sin(aoa/2), 0});
    SX q_ssa = SX::vertcat({cos(ssa/2), 0, 0, sin(-ssa/2)});
    SX q_BV = quatmul(q_aoa, q_ssa);

    // total drag coefficient
    SX CD = CD0_tot + pow(CL0 + CLa_tot * aoa,2) / (M_PI * e_o * AR);

    // Hydrodynamic Forces
    // Masterthesis equations. Assume Q = rho/2*S
    // Z  = (CL(aoa) + CLq * W(2) + CLaoa_dot(aoa_dot) + CLde*dE + CLdf*dF)*Q*V^2
    SX LIFT = (CL0 + CLa_tot * aoa + CZdf * dF) * dyn_press * S + (0.25 * CLq * c * S * rho) * V * W(1);

    // X  = (CD(aoa) + CXdf * dF) * Q * V^2
    SX DRAG = (CD + CXdf * dF) * dyn_press * S;

    // Y  = (CYb * ssa + CYp*W(1) + CYr* W(3) + CYdr * dR)*Q*V^2
    SX SF = (CYb * ssa + CYdr*dR) * dyn_press * S + 0.25 * (CYr * W(2) + CYp * W(0)) * (b * rho * S) * V;

    // Force return value
    Fhbrf = quatrot(q_BV, SX::vertcat({-DRAG, 0, -LIFT})) + SX::vertcat({0, SF, 0}); // why is SF already in brf?

    // Moments
    // Rolling Aerodynamic Moment
    // L = (Clb * ssa + Clp * W(1) + Clr * W(3) + Clda * dA + Cldr * dR) * Q * V^2 * b
    SX L = (Cl0 + Clb * ssa + Clda*dA + Cldr*dR) * dyn_press * S * b
           + (Clr * W(2) + Clp * W(0)) * (0.25 * rho * pow(b,2) * S) * V;

    // Pitching Aerodynamic Moment
    // M = (Cm0 + Cma*aoa + Cmaoa_dot*aoa_dot + Cmq*W(2) + Cmde*dE + Cmdf*dF) * Q * V^2* c
    SX M = (Cm0 + Cma * aoa  + Cmdf * dF) * dyn_press * S * c
           + Cmq * (0.25 * S * pow(c,2) * rho) * W(1) * V;

    // Yawing Aerodynamic Moment
    // N = (Cnb * ssa + Cnp * W(1) + Cnr * W(3) + Cnda * dA + Cndr * dR) * Q * V^2 * b
    SX N = (Cn0 + Cnb * ssa  + Cnda * dA + Cndr*dR) * dyn_press * S * b
           + (Cnp * W(0) + Cnr * W(2)) * (0.25 * S * pow(b,2) * rho) * V;

    // Angular motion equation in BRF
    SX Mhvrf = SX::vertcat({L, M, N});

    // Moment return value
    Mhbrf = quatrot(q_BV, Mhvrf);
}

void BoatDynamics::Propulsion(const SX &state, const SX &control, const BoatProperties &prop, SX &Ftbrf, SX &Mtbrf)
{
    // TODO: move constants into BoatProperties

    SX thrust = control(3);

    // propmode: "endurance"
    // double Ftbrfx = 71;   // [N] in the Steady State
    // double Mtbrfx = 6.18; // [Nm], optimal torque (qopt). Prop rotates cw wrt brf_x => reaction

    // propmode: "speed"
    double Ftbrfx = 93;   // [N] in the Steady State
    double Mtbrfx = 4.77; // [Nm], optimal torque (qopt). Prop rotates cw wrt brf_x => reaction
    Ftbrf = SX::vertcat({Ftbrfx, 0, 0});
    Ftbrf = thrust * Ftbrf;

    // Parasite/Reaction Torque: Moment around y due to thrust
    // d_prop_vec = [0.6766 0.7716 0.8666]
    // boat_sims.cbdepth = 1
    // double d_prop_cog = d_prop_vec(boat_sims.cbdepth);
    double d_prop_cog = 0.7716;
    double Mrbrf = d_prop_cog*Ftbrfx;

    Mtbrf = SX::vertcat({Mtbrfx, Mrbrf, 0});
    Mtbrf = thrust * Mtbrf;
}

void quat_error_mult(const SX &a, const SX &qr, SX &q, SX &eta)
{
    // Attitude error transfer to reference quaternion q
    const double f = 1;
    eta = 1 - SX::dot(a, a) / 8;
    SX nu = a / 2;
    SX dq = f * SX::vertcat({eta, nu});

    q = quatmul(dq, q);
    q /= SX::norm_2(q); // enforce unit norm constraint
}

BoatDynamics::EstimatorJacobian(const SX &state, const SX &control, const BoatProperties &prop)
{
    double t_samp = prop.estimator.t_samp;

    SX a = SX::sym("a", 3);      // Attitude error parametrisation
    SX bg = SX::sym("bg", 3);    // IMU accelerometer bias
    SX ba = SX::sym("ba", 3);    // IMU gyroscope bias

    SX u = SX::sym("u", 3);      // Control: Flaps, Ailerons, Rudder
    SX qr = SX::sym("qr", 4);     // reference, true of last iteration

    // Erroneous Quaternion
    SX q, eta;
    quat_error_mult(a, qr, q, eta);

    SX xe = SX::vertcat({v, W, r, a, bg, ba}); // estimator state
    SX xs = SX::vertcat({v, W, r, q});

    // (M)EKF dynamics/ Composition of total Jacobian

    // Model Jacobian to get linearized model
    // auto A_sys = jacobian_func(SXVector{xs, u});
    SX A_sys = jacobian;

    // df/da = df/dq*dq/da; Jq := dq/da
    SX Jq = SX::jacobian(q, a);

    // Matrix valued chain rule
    SX Aa = A_sys(Slice(0, 9), Slice(9, 13)) * Jq;

    // Attitude Dynamics (from Nokland)
    SX f_att_nl = SX::vertcat({
        0.5*(SX::eye(3) * eta + skew_mat(a)) * (W - bg),
        SX::zeros(3,1)
    });

    // Attitude Dynamics (from Markley)
    SX f_att_mk = SX::vertcat({
        -bg - skew_mat(W)*a,
        SX::zeros(3,1)
    });

    SX Jatt = SX::jacobian(f_att_mk, xe);

    SX A = SX::vertcat({
        SX::horzcat({A_sys(Slice(0, 9), Slice(0, 9)), Aa, SX::zeros(9, 6)}),
        Jatt,
        SX::zeros(3, xe.size1())
    });
    // A_func = Function("A_func", {xe, u, qr}, {F});

    // Discretization
    SX F = SX::eye(xe.size1()) + t_samp*A; // Euler
    Function F_func = Function("F_func", {xe, u, qr}, {F});

    // this->Estimator.F = F;
    // this->Estimator.F_func = F_func;

    // Outputmap and Sensitivity
    SX g = SX::vertcat({0, 0, -prop.env.g}); // in BRF: z = up
    SX r_ant = SX::vertcat({
        prop.sensor.r_ant[0],
        prop.sensor.r_ant[1],
        prop.sensor.r_ant[2]
    });
    SX g_BRF = quatrot(q, g);

    // Velocity is in BRF
    SX y_vgns = quatrot_inverse(q, v) + quatrot_inverse(q, skew_mat(W) * r_ant);
    SX y_wgyro = W + bg;
    SX y_rgns = r + quatrot_inverse(q, r_ant);

    //y_acc = v_dot_BRF - quatrot(q,[0; g]) + ba;
    SX y_acc = -quatrot(q, SX::vertcat({0,g})) + ba;

    // Complete Output-map h(x) and Output Sensitivity H
    SX h = SX::vertcat({y_vgns, y_wgyro, y_rgns, y_acc});
    Function h_func = Function("h_func", {xe,u, qr, r_ant}, {h});
    SX H = SX::jacobian(h, xe);
    Function H_func = Function("H_func", {xe,u, qr, r_ant}, {H});

    // this->Estimator.h = h;
    // this->Estimator.h_func = h_func;
    // this->Estimator.H = H;
    // this->Estimator.H_func = H_func;
}

BoatDynamics::BoatDynamics(const BoatProperties &prop)
{
    double g = prop.env.g;
    double mass = prop.inertia.mass;
    double mass_cad = prop.inertia.mass_cad;
    double mcorr = mass / mass_cad; // Until better values: make boat uniformly heavier/inert

    double Ixy = mcorr * prop.inertia.Ixy;
    double Ixz = mcorr * prop.inertia.Ixz;
    double Iyz = mcorr * prop.inertia.Iyz;

    double Ixx = mcorr * prop.inertia.Ixx;
    double Iyy = mcorr * prop.inertia.Iyy;
    double Izz = mcorr * prop.inertia.Izz;

    SX Jbrf = SX::vertcat({  // [kg m^2], inertial tensor
        SX::horzcat({Ixx, Ixy, Ixz}),
        SX::horzcat({Ixy, Iyy, Iyz}),
        SX::horzcat({Ixz, Iyz, Izz})
    });

    // State and Control variables
    SX v    = SX::sym("v", 3);   // [m/s] translational velocity of the CoG in BRF
    SX W    = SX::sym("W", 3);   // [rad/s] angular velocities of boat in BRF
    SX r    = SX::sym("r", 3);   // [m] position of the CoG in IRF
    SX q_BI = SX::sym("q", 4);   // unit quaternion for transformation from IRF to BRF

    SX dF   = SX::sym("dF");     // Flaps
    SX dA   = SX::sym("dA");     // Aileron deflection [reserved, but not used] [rad]
    SX dR   = SX::sym("dR");     // Rudder deflection [rad]
    // SX dE   = SX::sym("dE");     // Elevator deflection [positive down] [rad]
    SX thrust = SX::sym("thrust"); // Thrust [-] (between 0.0 and 1.0)

    SX state    = SX::vertcat({v, W, r, q_BI});
    SX control  = SX::vertcat({dF, dA, dR, thrust}); // TODO: why is elevator deflection dE not used?

    // Gravity
    SX Fgbrf = mass*quatrot(q_BI, SX::vertcat({0,0,g})); // inertial reference frame, NED

    // Propulsion
    SX Ftbrf, Mtbrf;
    Propulsion(state, control, prop, Ftbrf, Mtbrf);

    // Hydrodynamics
    SX Fhbrf, Mhbrf, aoa, ssa; // TODO: aoa, ssa unused?
    Hydrodynamics(state, control, prop, Fhbrf, Mhbrf, aoa, ssa);

    // Buoyancy
    // SX Fbbrf, Mbbrf;
    // TODO

    // Damping
    SX Df = SX::diag(SX::vertcat({-10, -20, -5})); // TODO: parametrize
    SX Fdbrf = SX::mtimes(Df, v);

    SX Dm = SX::diag(SX::vertcat({-2000, -800, -300}));
    SX Mdbrf = SX::mtimes(Dm, W);


    SX Fbrf = Fhbrf + Ftbrf + Fdbrf + Fgbrf; // + Fbbrf
    SX Mbrf = Mhbrf + Mdbrf; // + Mtbrf + Mbbrf; // TODO: why no thrust moment?

    // Boat translational velocity in BRF
    SX v_dot_brf = Fbrf / mass - SX::cross(W,v);

    // Boat roataional velocity in BRF
    SX W_dot_brf = SX::mtimes(SX::inv(Jbrf), (Mbrf - SX::cross(W, SX::mtimes(Jbrf, W))));

    // Dynamic Equations: Kinematics
    SX r_dot_irf = quatrot_inverse(q_BI, v);
    SX q_BI_dot  = 0.5 * quatmul(q_BI, SX::vertcat({0,W})); // TODO: verify correctness

    // Differential equation
    SX dynamics = SX::vertcat({v_dot_brf, W_dot_brf, r_dot_irf, q_BI_dot});

    Function dynamics_func = Function("dynamics", {state, control}, {dynamics});
    SX jacobian = SX::jacobian(dynamics, state);
    Function jacobian_func = Function("jacobian", {state, control}, {jacobian});

    // define RK4 integrator scheme
    SX X = SX::sym("X", state.size1());
    SX U = SX::sym("U", control.size1());
    SX dT = SX::sym("dT");

    // get symbolic expression for RK4 integrator
    SX integrator = rk4_symbolic(X, U, dynamics_func, dT);
    Function integrator_func = Function("RK4", {X,U,dT},{integrator});

    // assign class attributes
    this->State = state;
    this->Control = control;

    this->SymDynamics = dynamics;
    this->SymIntegartor = integrator;
    this->SymJacobian = jacobian;

    this->NumDynamics = dynamics_func;
    this->NumJacobian = jacobian_func;
    this->NumIntegrator = integrator_func;

    EstimatorJacobian(state, control, prop);
}

} // namespace bifoiler
