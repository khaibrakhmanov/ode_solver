#include "../lib/ode_solver.h"
#include "eqs.h"

/// @brief class: Runge-Kutta solver for a single ODE
class TRungeKuttaSolver : public TBaseRungeKuttaSolver
{
private:
    /// @brief ODE to be solved (to be initialised)
    TSingleODE *ode;

public:
    TRungeKuttaSolver(/* args */);
    ~TRungeKuttaSolver();

    /// @brief Import initialised ODE into the solver
    /// @param ode - single ODE to be solved
    void ImportODE(TSingleODE *ode);

    /// @brief Calculate the RK-coefficients
    /// @param y0 - vector of the initial values of the independent varables
    /// @param dt - time step
    /// @param t - time
    void CalcCoefficients(const TArray y0, const double dt, const double t);

    void UserCheck() {};
};


/// @brief class: Runge-Kutta solver for a single ODE
class TRungeKuttaSolverSystem : public TBaseRungeKuttaSolver
{
private:
    /// @brief ODEs to be solved (to be initialised)
    TODEs* odes;

public:
    TRungeKuttaSolverSystem(/* args */);
    ~TRungeKuttaSolverSystem();

    /// @brief Import initialised ODE into the solver
    /// @param ode - single ODE to be solved
    void ImportODEs(TODEs* ode);

    /// @brief Calculate the RK-coefficients
    /// @param y0 - vector of the inital values of the independent varable
    /// @param dt - time step
    /// @param t - time
    void CalcCoefficients(const TArray y0, const double dt, const double t);

    void UserCheck() {};
};
