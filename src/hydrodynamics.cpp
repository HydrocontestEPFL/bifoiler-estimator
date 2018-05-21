#include <iostream>
#include <math.h>

#include <casadi/casadi.hpp>
#include <yaml-cpp/yaml.h>

using namespace casadi;

namespace bifoiler {

struct BoatProperties
{
    std::string name;

    struct {
        double rho_sh2o;
    } env;

    struct {
        double mass;
        double mas_cad;
        double Ibxx_xflr5;
        double Ibyy_xflr5;
        double Ibzz_xflr5;
        double Ibxz_xflr5;
        double Ixy;
        double Ixz;
        double Iyz;
        double Ixx;
        double Iyy;
        double Izz;
    } inertia;

    struct {
        double ARff;
        double areaff;
        double mac;
        double wingspanff;
    } foils;

    struct {
        double CL0;
        double CLa_total;
        double e_oswald;
        double CD0_total;
        double CYb;
        double Cm0;
        double Cma;
        double Cn0;
        double Cnb;
        double Cl0;
        double Clb;
        double CLq;
        double Cmq;
        double CYr;
        double Cnr;
        double Clr;
        double CYp;
        double Clp;
        double Cnp;
        double CXdf;
        double CYdr;
        double CZde;
        double CZdf;
        double CLda;
        double CLdr;
        double CMde;
        double CMdf;
        double CNda;
        double CNdr;
    } hydrodynamic;
};

class BoatDynamics
{
public:
    casadi::SX State;
    casadi::SX Control;
    casadi::SX SymDynamics;
    casadi::SX SymIntegartor;
    casadi::SX SymJacobian;

    casadi::Function NumDynamics;
    casadi::Function NumIntegrator;
    casadi::Function NumJacobian;

    BoatDynamics(const BoatProperties &boat_prop);
    static void Hydrodynamics(SX &state, SX &control, const BoatProperties &prop, SX &Fhbrf, SX &Mhbrf, SX &aoa, SX &ssa);
    static void Propulsion(const BoatProperties &prop, SX &Ftbrf, SX &Mtbrf)
    static BoatProperties LoadProperties(const std::string &filename);
};

BoatProperties BoatDynamics::LoadProperties(const std::string &filename)
{
    YAML::Node config = YAML::LoadFile(filename);

    BoatProperties prop;
    prop.name = config["name"].as<std::string>();
    prop.env.rho_sh2o = config["env"]["rho_sh2o"].as<double>();

    auto foils = config["foils"];
    prop.foils.ARff = foils["ARff"].as<double>();
    prop.foils.areaff = foils["areaff"].as<double>();
    prop.foils.mac = foils["mac"].as<double>();
    prop.foils.wingspanff = foils["wingspanff"].as<double>();

    auto hydro = config["hydrodynamic"];
    prop.hydrodynamic.CL0 = hydro["CL0"].as<double>();
    prop.hydrodynamic.CLa_total = hydro["CLa_total"].as<double>();
    prop.hydrodynamic.e_oswald = hydro["e_oswald"].as<double>();
    prop.hydrodynamic.CD0_total = hydro["CD0_total"].as<double>();
    prop.hydrodynamic.CYb = hydro["CYb"].as<double>();
    prop.hydrodynamic.Cm0 = hydro["Cm0"].as<double>();
    prop.hydrodynamic.Cma = hydro["Cma"].as<double>();
    prop.hydrodynamic.Cn0 = hydro["Cn0"].as<double>();
    prop.hydrodynamic.Cnb = hydro["Cnb"].as<double>();
    prop.hydrodynamic.Cl0 = hydro["Cl0"].as<double>();
    prop.hydrodynamic.Clb = hydro["Clb"].as<double>();
    prop.hydrodynamic.CLq = hydro["CLq"].as<double>();
    prop.hydrodynamic.Cmq = hydro["Cmq"].as<double>();
    prop.hydrodynamic.CYr = hydro["CYr"].as<double>();
    prop.hydrodynamic.Cnr = hydro["Cnr"].as<double>();
    prop.hydrodynamic.Clr = hydro["Clr"].as<double>();
    prop.hydrodynamic.CYp = hydro["CYp"].as<double>();
    prop.hydrodynamic.Clp = hydro["Clp"].as<double>();
    prop.hydrodynamic.Cnp = hydro["Cnp"].as<double>();
    prop.hydrodynamic.CXdf = hydro["CXdf"].as<double>();
    prop.hydrodynamic.CYdr = hydro["CYdr"].as<double>();
    prop.hydrodynamic.CZde = hydro["CZde"].as<double>();
    prop.hydrodynamic.CZdf = hydro["CZdf"].as<double>();
    prop.hydrodynamic.CLda = hydro["CLda"].as<double>();
    prop.hydrodynamic.CLdr = hydro["CLdr"].as<double>();
    prop.hydrodynamic.CMde = hydro["CMde"].as<double>();
    prop.hydrodynamic.CMdf = hydro["CMdf"].as<double>();
    prop.hydrodynamic.CNda = hydro["CNda"].as<double>();
    prop.hydrodynamic.CNdr = hydro["CNdr"].as<double>();

    return prop;
}

SX quatmul(const SX &q1, const SX &q2)
{
    SX s1 = q1(0);
    SX v1 = q1(Slice(1,4));

    SX s2 = q2(0);
    SX v2 = q2(Slice(1,4));

    SX q = SX::vertcat({s1*s2 - SX::dot(v1,v2), SX::cross(v1,v2) + s1 * v2 + s2 * v1});
    return q;
}

SX quatinv(const SX &q)
{
    SX inv_q = SX::vertcat({q(0), -q(1), -q(2), -q(3)}) / SX::norm_2(q);
    return inv_q;
}

SX quatrot(const SX &q, const SX &r)
{
    SX qinv = quatinv(q);
    SX qr   = SX::vertcat({0, r});
    SX qrr  = quatmul(quatmul(qinv, qr),q);
    return qrr(Slice(1,4));
}

void BoatDynamics::Hydrodynamics(SX &state, SX &control, const BoatProperties &prop, SX &Fhbrf, SX &Mhbrf, SX &aoa, SX &ssa)
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
    SX dE = control(3); // Elevator

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

void BoatDynamics::Propulsion(const BoatProperties &prop, SX &Ftbrf, SX &Mtbrf)
{
    // TODO: move constants into BoatProperties

    // propmode: "endurance"
    // double Ftbrfx = 71;   // [N] in the Steady State
    // double Mtbrfx = 6.18; // [Nm], optimal torque (qopt). Prop rotates cw wrt brf_x => reaction

    // propmode: "speed"
    double Ftbrfx = 93;   // [N] in the Steady State
    double Mtbrfx = 4.77; // [Nm], optimal torque (qopt). Prop rotates cw wrt brf_x => reaction
    Ftbrf = SX::vertcat({Ftbrfx, 0, 0});

    // Parasite/Reaction Torque: Moment around y due to thrust
    // d_prop_vec = [0.6766 0.7716 0.8666]
    // boat_sims.cbdepth = 1
    // double d_prop_cog = d_prop_vec(boat_sims.cbdepth);
    double d_prop_cog = 0.7716;
    double Mrbrf = d_prop_cog*Ftbrfx;

    Mtbrf = SX::vertcat({Mtbrfx, Mrbrf, 0});
}

BoatDynamics::BoatDynamics(const BoatProperties &prop)
{
#if 1
    double g = boat_pm.env.g;
    double mtot = boat_pm.inertia.mass;

    // TODO
    // Jbrf = [boat_pm.Ixx boat_pm.Ixy boat_pm.Ixz; boat_pm.Ixy boat_pm.Iyy boat_pm.Iyz; boat_pm.Ixz boat_pm.Iyz boat_pm.Izz]; %[kg m^2], inertial tensor
    // Jbrf = Jbrf./boat_pm.mtot_cad*mtot;              % Until better values: make boat uniformly heavier/inert

    // State and Control variables
    SX v    = SX.sym("v", 3);   // [m/s] translational velocity of the CoG in BRF
    SX W    = SX.sym("W", 3);   // [rad/s] angular velocities of boat in BRF
    SX r    = SX.sym("r", 3);   // [m] position of the CoG in IRF
    SX q_BI = SX.sym("q", 4);   // unit quaternion for transformation from IRF to BRF

    SX dF   = SX.sym("dF");     // Flaps
    SX dA   = SX.sym("dA");     // Aileron deflection [reserved, but not used] [rad]
    SX dR   = SX.sym("dR");     // Rudder deflection [rad]
    SX dE   = SX.sym("dE");     // Elevator deflection [positive down] [rad]

    SX state    = SX::vertcat({v, W, r, q_BI});
    SX control  = SX::vertcat({dF, dA, dR, dE});

    // Gravity
    SX Fgbrf = quatrot(q_BI, SX::vertcat({0,0,g}))

    // Propulsion
    SX Ftbrf, Mtbrf;
    Propulsion(prop, Ftbrf, Mtbrf);

    // Hydrodynamics
    SX Fhbrf, Mhbrf;
    Hydrodynamics(prop, state, control, Fhbrf, Mhbrf);

    // Buoyancy
    // SX Fbbrf, Mbbrf;
    // TODO

    // Damping
    Df = SX::diag({-10, -20, -5}); // TODO: parametrize
    Fdbrf = Df*v;

    Dm = SX::diag({-2000, -800, -300});
    Mdbrf = Dm*W;


    SX Fbrf = Fhbrf + Ftbrf + Fdbrf + Fgbrf; // + Fbbrf
    SX Mbrf = Mhbrf + Mtbrf + Mdbrf; // + ...

    // Boat translational velocity in BRF
    v_dot_brf = Fbrf/mtot - cross(W,v);

    // Boat roataional velocity in BRF
    W_dot_brf = inv(Jbrf) * (Mbrf - cross(W, Jbrf * W)); // TODO; Jbrf

    // Dynamic Equations: Kinematics
    qv        = quatmul(q_BI, SX::vertcat({0; v}));
    qv_q      = quatmul(qv, quatinv(q_BI));
    r_dot_irf = qv_q(Slice(1:4));
    q_BI_dot  = 0.5 * quatmul(q_BI, SX::vertcat({0,W}));

    // Differential equation
    SX dynamics = SX::vertcat({v_dot_brf, W_dot_brf, r_dot_irf, q_BI_dot});

    dynamics_func = Function("dynamics", {state,control}, {dynamics});
    jacobian = dynamics_func.jacobian(state);
    jacobian_func = Function('jacobian', {state, control}, {jacobian});

    // define RK4 integrator scheme
    SX X = SX::sym("X", 13);
    SX U = SX::sym("U", 3);
    SX dT = SX::sym("dT");

    // get symbolic expression for RK4 integrator
    SX integrator = kmath::rk4_symbolic(X, U, dynamics_func, dT);
    Function integrator_func = Function("RK4", {X,U,dT},{sym_integrator});

    // assign class attributes
    this->State = state;
    this->Control = control;

    this->SymDynamics = dynamics;
    this->SymIntegartor = integrator;

    this->NumDynamics = dynamics_func;
    this->NumJacobian = jacobian_func;
    this->NumIntegrator = integrator_func;
#endif

    double t_samp = boat_pm.estimator.t_samp; // TODO

    SX a = SX::sym('a',3);      // Attitude error parametrisation
    SX bg = SX::sym('bg',3);    // IMU accelerometer bias
    SX ba = SX::sym('ba',3);    // IMU gyroscope bias

    // Attitude error transfer to reference quaternion q
    const double f = 1;
    SX eta = 1 - SX::dot(a,a) / 8;
    SX nu = a / 2;
    SX dq = f * SX::vertcat({eta, nu})

    SX q = quatmul(dq, q);
    q /= q.norm_2(); // enforce unit norm constraint

    xe = SX::vertcat({v, W, r, a, bg, ba}); // estimator state
    xs = SX::vertcat({v, W, r, q});

    // (M)EKF dynamics/ Composition of total Jacobian

    A = SX::vertcat({
        SX::horzcat({A_sys(Slice(0,9), Slice(0,9)), Aa, SX::zeros(9,6)}),
        Jatt,
        SX::zeros(3, xe.size1())
    });

    // Discretization
    F = SX::eye(xe.size1()) + t_samp*A; // Euler
    F_func = Function("F_func", {xe, u, qr}, {F});
}

} // namespace bifoiler

using namespace bifoiler;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " path/to/config.yaml" << std::endl;
        return -1;
    }

    std::string config_file(argv[1]);
    auto prop = BoatDynamics::LoadProperties(config_file);

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
