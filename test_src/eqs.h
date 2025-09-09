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

    double RHS(const double t, const double y) const;
    double Solution(const double t) const;
};

/// @brief class: ODE with the linear derivative: dy/dt = a - b*y
class TLinODE : public TSingleODE
{
public:
    /// @brief 'Exponential' ODE constructor
    /// @param a
    /// @param b
    /// @param ic
    TLinODE(const double a, const double b, const double ic);
    ~TLinODE();

    double RHS(const double t, const double y) const;
    double Solution(const double t) const;
};

/// @brief class: ODE with the quadratic derivative: dy/dt = a - b*y^2
class TQuadODE : public TSingleODE
{
public:
    /// @brief 'Quadratic' ODE constructor
    /// @param a
    /// @param b
    /// @param ic
    TQuadODE(const double a, const double b, const double ic);
    ~TQuadODE();

    double RHS(const double t, const double y) const;
    double Solution(const double t) const;
};