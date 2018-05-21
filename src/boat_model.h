#ifndef BOAT_MODEL_H
#define BOAT_MODEL_H

#include <casadi/casadi.hpp>
#include "boat_properties.h"

namespace bifoiler {

class BoatDynamics
{
public:
    BoatDynamics(const BoatProperties &boat_prop);
    static void Hydrodynamics(const casadi::SX &state,
                              const casadi::SX &control,
                              const BoatProperties &prop,
                              casadi::SX &Fhbrf,
                              casadi::SX &Mhbrf,
                              casadi::SX &aoa,
                              casadi::SX &ssa);
    static void Propulsion(const casadi::SX &state,
                           const casadi::SX &control,
                           const BoatProperties &prop,
                           casadi::SX &Ftbrf,
                           casadi::SX &Mtbrf);

private:
    casadi::SX State;
    casadi::SX Control;
    casadi::SX SymDynamics;
    casadi::SX SymIntegartor;
    casadi::SX SymJacobian;

    casadi::Function NumDynamics;
    casadi::Function NumIntegrator;
    casadi::Function NumJacobian;
};

} // namespace bifoiler

#endif /* BOAT_MODEL_H */
