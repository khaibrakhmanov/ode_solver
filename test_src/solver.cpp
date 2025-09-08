#include "solver.h"

TRungeKuttaSolver::TRungeKuttaSolver(/* args */)
{
}

TRungeKuttaSolver::~TRungeKuttaSolver()
{
}

void TRungeKuttaSolver::ImportODE(TSingleODE *ode)
{
    this->ode = ode;
}

void TRungeKuttaSolver::CalcCoefficients(const TArray v0, const double dt, const double t)
{
    TArray2D* coeffs = GetCoefficientsPtr();

    (*coeffs)[0][0] = ode->RHS(t, v0[0]);
    (*coeffs)[1][0] = ode->RHS(t + 0.5 * dt, v0[0] + 0.5 * dt * GetC(0, 0));
    (*coeffs)[2][0] = ode->RHS(t + 0.5 * dt, v0[0] + 0.5 * dt * GetC(1, 0));
    (*coeffs)[3][0] = ode->RHS(t + dt, v0[0] + dt * GetC(2, 0));
}
