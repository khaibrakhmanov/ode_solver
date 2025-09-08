#pragma once
#include "ode.h"


/// @brief class: ODE with the exponentially decaying solution: dy/dt = -alpha*y
class TExpODE : public TSingleODE
{
public:
    /// @brief 'Exponential' ODE constructor
    /// @param alpha decay rate (stiffness parameter)
    /// @param ic initial condition
    TExpODE(const double alpha, const double ic);
    ~TExpODE();

    double RHS(const double t, const double x) const;
    double Solution(const double t, const double x) const;
};
