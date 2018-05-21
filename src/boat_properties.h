#ifndef BOAT_PROPERTIES_H
#define BOAT_PROPERTIES_H

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

    static BoatProperties Load(const std::string &filename);
    // static void Save(const std::string &filename); // TODO
};

#endif /* BOAT_PROPERTIES_H */
