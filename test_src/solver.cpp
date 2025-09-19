#include "solver.h"

TRungeKuttaSolverSystem::TRungeKuttaSolverSystem(/* args */)
{
}

TRungeKuttaSolverSystem::~TRungeKuttaSolverSystem()
{
}

void TRungeKuttaSolverSystem::ImportODEs(TODEs* ode)
{
    this->odes = ode;
}

void TRungeKuttaSolverSystem::CalcCoefficients(const TArray y0, const double dt, const double t)
{
    TArray2D* coeffs = GetCoefficientsPtr();

    auto N = y0.size();
    auto rhs = odes->GetRHSPtr();

    y1.clear();
    for (size_t i = 0; i < N; i++)
    {
        (*coeffs)[0][i] = rhs->at(i)(t, y0);
        y1.push_back(y0[i] + 0.5 * dt * GetC(0, i));
    }

    y2.clear();
    for (size_t i = 0; i < N; i++)
    {
        (*coeffs)[1][i] = rhs->at(i)(t + 0.5 * dt, y1);
        y2.push_back(y0[i] + 0.5 * dt * GetC(1, i));
    }

    y3.clear();
    for (size_t i = 0; i < N; i++)
    {
        (*coeffs)[2][i] = rhs->at(i)(t + 0.5 * dt, y2);
        y3.push_back(y0[i] + dt * GetC(2, i));
    }

    for (size_t i = 0; i < N; i++)
        (*coeffs)[3][i] = rhs->at(i)(t + dt, y3);
}
