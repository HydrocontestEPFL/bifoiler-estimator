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

BoatProperties LoadProperties(const std::string &filename)
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

void BoatHydrodynamics(const BoatProperties &prop)
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
    SX v = SX::sym("v", 3); // Translational Velocity/Attitude towards Velocity Reference Frame
    // SX v = SX::vertcat(state(Slice(0,3)));
    SX W = SX::sym("W", 3); // Rotational Rate/Attitude Derivative
    // SX W = state(Slice(3,6));

    // Control variables definition
    SX dF = SX::sym("dF"); // Flaps
    SX dA = SX::sym("dA"); // Aileron
    SX dR = SX::sym("dR"); // Rudder
    SX dE = SX::sym("dE"); // Elevator

    // Variables and Orientation
    const double v1_eps   = 0.0001;    // Enable v(1)_0 = 0 w/o getting NaN
    SX V = SX::norm_2(v)+ v1_eps;        // Absolute Value of Velocity

    SX ssa = asin(v(1) / V);               // side slip angle [rad]
    SX aoa = atan2(v(2) , v(0) + v1_eps);  // angle of attack definition [rad]

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

    SX Fhbrf = quatrot(q_BV, SX::vertcat({-DRAG, 0, -LIFT})) + SX::vertcat({0, SF, 0}); // why is SF already in brf?

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
    SX Mhbrf = quatrot(q_BV, Mhvrf);
    
    SX state    = SX::vertcat({v,w});
    SX control  = SX::vertcat({...});
    SX dynamics = SX::vertcat({Fhbrf,Mhbrf});
    
    Function dynamics = Function("dynamics",{state,control},{dynamics});
    dynamics.generate();

    std::cout << Fhbrf << std::endl;
    std::cout << Mhbrf << std::endl;
}

} // namespace bifoiler

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " path/to/config.yaml" << std::endl;
        return -1;
    }

    std::string config_file(argv[1]);
    auto prop = bifoiler::LoadProperties(config_file);
    bifoiler::BoatHydrodynamics(prop);
}
