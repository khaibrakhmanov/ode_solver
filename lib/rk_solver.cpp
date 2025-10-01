#include "rk_solver.h"

TRungeKuttaSolver::TRungeKuttaSolver()
{
	richardsonExtrapolation.TurnOff();
}

void TRungeKuttaSolver::ImportODEs(TODEs *ode)
{
	this->odes = ode;
}

void TRungeKuttaSolver::CalcCoefficients(const TArray y0, const double dt, const double t)
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

double TRungeKuttaSolver::GetC(const size_t i, const size_t j)
{
	return coeff[i][j];
};

void TRungeKuttaSolver::SetC(const size_t i, const size_t j, const double arg)
{
	coeff[i][j] = arg;
};

TArray2D* TRungeKuttaSolver::GetCoefficientsPtr()
{
	return &coeff;
}

void TRungeKuttaSolver::AllocateCoefficients()
{
	AllocateVector2D(&coeff, GetOrder(), GetVariablesNumber());
};

size_t TRungeKuttaSolver::GetOrder()
{
	return order;
};

void TRungeKuttaSolver::SetOrder(const size_t arg)
{
	order = arg;
};

void TRungeKuttaSolver::CalcIntemediateVariables(const TArray y0)
{
    for(std::size_t i = 0; i < GetVariablesNumber(); i++)
    {
        y1[i] = y0[i] + 0.5 * dt * GetC(0, i);
        y2[i] = y0[i] + 0.5 * dt * GetC(1, i);
        y3[i] = y0[i] + dt * GetC(2, i);
    }
}

void TRungeKuttaSolver::Step()
{
		CalcCoefficients(initialVector, GetTimeStep(), GetTime());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(sol, initialVector)
	{
#pragma omp for nowait
		for (int i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			sol[i_sol] = initialVector[i_sol] + GetTimeStep() / 6.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() / 6.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif
}

void TRungeKuttaSolver::Step2()
{
		CalcCoefficients(initialVector, 0.5 * GetTimeStep(), GetTime());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(solHalf, initialVector)
	{
#pragma omp for nowait
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			solHalf[i_sol] = initialVector[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		solHalf[i_sol] = initialVector[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

	CalcCoefficients(solHalf, 0.5 * GetTimeStep(), GetTime() + 0.5 * GetTimeStep());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(sol2, solHalf)
	{
#pragma omp for nowait
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			sol2[i_sol] = solHalf[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol2[i_sol] = solHalf[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif
}

//-------- Runge-Kutta-Fehlberg 4(5) solver -----------------------
TBaseRungeKuttaFehlberg45Solver::TBaseRungeKuttaFehlberg45Solver()
{
	SetOrder(5);
	b1_4 = 25.0 / 216.0;
	b3_4 = 1408.0 / 2565.0;
	b4_4 = 2197.0 / 4164.0;
	b5_4 = -0.2;

	b1_5 = 16.0 / 135.0;
	b3_5 = 6656.0 / 12825.0;
	b4_5 = 28561.0 / 56430.0;
	b5_5 = -9.0 / 50.0;
	b6_5 = 2.0 / 55.0;
}

void TBaseRungeKuttaFehlberg45Solver::Step()
{
	CalcCoefficients(initialVector, GetTimeStep(), GetTime());

	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() * (b1_4 * GetC(0, i_sol) + b3_4 * GetC(2, i_sol) + b4_4 * GetC(3, i_sol) - b5_4 * GetC(4, i_sol));
	};
}

void TBaseRungeKuttaFehlberg45Solver::Step2()
{
	CalcCoefficients(initialVector, GetTimeStep(), GetTime());

	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol2[i_sol] = initialVector[i_sol] + GetTimeStep() * (b1_5 * GetC(0, i_sol) + b3_5 * GetC(2, i_sol) + b4_5 * GetC(3, i_sol) + b5_5 * GetC(4, i_sol) + b6_5 * GetC(5, i_sol));
	};
}