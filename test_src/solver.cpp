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

// ----------------------------------------------------------------------------------------

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

void TRungeKuttaSolverSystem::CalcCoefficients(const TArray v0, const double dt, const double t)
{
    TArray2D* coeffs = GetCoefficientsPtr();

    auto N = v0.size();
    auto rhs = odes->GetRHSPtr();

    for (size_t i = 0; i < N; i++)
    {
        (*coeffs)[0][i] = (rhs->at(i))(t, v0);
        // сумма вектора и функции?... сделать C(0, i) векторной функцей7?
        (*coeffs)[1][i] = (rhs->at(i))(t + 0.5 * dt, v0[i] + 0.5 * dt * GetC(0, i));

        (*coeffs)[0][i] = ode->RHS(t, v0[i]);
        (*coeffs)[1][i] = ode->RHS(t + 0.5 * dt, v0[i] + 0.5 * dt * GetC(0, i));
        (*coeffs)[2][i] = ode->RHS(t + 0.5 * dt, v0[i] + 0.5 * dt * GetC(1, i));
        (*coeffs)[3][i] = ode->RHS(t + dt, v0[i] + dt * GetC(2, i));
    }

}
