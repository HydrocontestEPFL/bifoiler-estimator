#ifndef ESTIMATOR_MODEL_H
#define ESTIMATOR_MODEL_H

#include <casadi/casadi.hpp>

#include "boat_model.h"
#include "boat_properties.h"

namespace bifoiler {

class EstimatorModel
{
public:
    EstimatorModel(BoatDynamics &dynamics, const BoatProperties &prop);

    casadi::SX F; // propagation matrix
    casadi::SX A; // Estimator jacobian
    casadi::SX h;
    casadi::SX H;
    casadi::Function F_func;
    casadi::Function A_func;
    casadi::Function h_func;
    casadi::Function H_func;
};

} // namespace bifoiler

#endif /* ESTIMATOR_MODEL_H */
