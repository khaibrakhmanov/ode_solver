#include "gear_solver.h"
#include <algorithm>
#include <cmath>

void TBaseGearSolver::SetOrder(const size_t _order)
{
	order = _order;

	switch (GetOrder())
	{
	case 2:
		c_u0 = 1.0 / 3.0;
		c_u1 = 4.0 / 3.0;
		c_f = 2.0 / 3.0;
		break;
	case 4:
		c_u0 = 3.0 / 25.0;
		c_u1 = 16.0 / 25.0;
		c_u2 = 36.0 / 25.0;
		c_u3 = 48.0 / 25.0;
		c_f = 12.0 / 25.0;
		break;
	case 6:
		c_u0 = 10.0 / 147.0;
		c_u1 = 72.0 / 147.0;
		c_u2 = 225.0 / 147.0;
		c_u3 = 400.0 / 147.0;
		c_u4 = 450.0 / 147.0;
		c_u5 = 360.0 / 147.0;
		c_f = 60.0 / 147.0;
		break;
	default:
		c_u0 = 3.0 / 25.0;
		c_u1 = 16.0 / 25.0;
		c_u2 = 36.0 / 25.0;
		c_u3 = 48.0 / 25.0;
		c_f = 12.0 / 25.0;
		break;
	}
};

size_t TBaseGearSolver::GetOrder()
{
	return order;
};

TBaseGearSolver::TBaseGearSolver()
{
	maxIterationsNumber = 20;
};

void TBaseGearSolver::AllocateCoefficients()
{
	AllocateVector1D(&u0, GetVariablesNumber());
	AllocateVector1D(&u1, GetVariablesNumber());
	AllocateVector1D(&u2, GetVariablesNumber());
	AllocateVector1D(&u3, GetVariablesNumber());
	AllocateVector1D(&u4, GetVariablesNumber());
	AllocateVector1D(&u5, GetVariablesNumber());
	AllocateVector1D(&u_s, GetVariablesNumber());
	AllocateVector1D(&u_s1, GetVariablesNumber());
	AllocateVector1D(&u_guess, GetVariablesNumber());
};

double TBaseGearSolver::GetExplicitPart(const size_t i_sol)
{
	switch (GetOrder())
	{
	case 2:
		return c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	case 4:
		return c_u3 * u3[i_sol] - c_u2 * u2[i_sol] + c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	case 6:
		return c_u5 * u5[i_sol] - c_u4 * u4[i_sol] + c_u3 * u3[i_sol] - c_u2 * u2[i_sol] + c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	default:
		return c_u3 * u3[i_sol] - c_u2 * u2[i_sol] + c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	}
};

double TBaseGearSolver::GetPhi(const size_t i_sol, const TArray u)
{
	return GetExplicitPart(i_sol) + c_f * GetTimeStep() * rhs(i_sol, GetTime(), u);
};

void TBaseGearSolver::ReadStartingPoints(const size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3)
{
	u0[i_sol] = _u0;
	u1[i_sol] = _u1;
	u2[i_sol] = _u2;
	u3[i_sol] = _u3;
};

void TBaseGearSolver::ReadStartingPoints(const size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3, const double _u4, const double _u5)
{
	u0[i_sol] = _u0;
	u1[i_sol] = _u1;
	u2[i_sol] = _u2;
	u3[i_sol] = _u3;
	u4[i_sol] = _u4;
	u5[i_sol] = _u5;
};

void TBaseGearSolver::ReadStartingPoints(const size_t i_sol, const double _u0, const double _u1)
{
	u0[i_sol] = _u0;
	u1[i_sol] = _u1;
};

void TBaseGearSolver::ReadStartingPoints(const TArray2D u_in)
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		for (size_t i = 0; i < GetOrder(); i++)
		{
			u_start[i_sol][i] = u_in[i_sol][i];
		};
	};
};

void TBaseGearSolver::ReadGuess(const size_t i_sol, const double u_g)
{
	u_guess[i_sol] = u_g;
}

void TBaseGearSolver::Step()
{
	TArray::iterator vmax_Iter;
	SetErrorMax(1.0);

	//initial iteration
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		switch (GetOrder())
		{
		case 2:
			u_s[i_sol] = u1[i_sol];
			break;
		case 4:
			if (u2[i_sol] > 1e-20)
				u_s[i_sol] = u3[i_sol] * u3[i_sol] / u2[i_sol];
			else
				u_s[i_sol] = 2.0 * u3[i_sol] - u2[i_sol];
			break;
		case 6:
			u_s[i_sol] = u5[i_sol] * u5[i_sol] / u4[i_sol];
			break;
		default:
			u_s[i_sol] = u3[i_sol] * u3[i_sol] / u2[i_sol];
			break;
		};
	};

	reachedIteration = 0;
	//iteration cycle
	while ((GetErrorMax() >= GetTolerance()) && (reachedIteration < maxIterationsNumber))
	{
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			u_s1[i_sol] = GetPhi(i_sol, u_s);
			errorVector[i_sol] = fabs(u_s1[i_sol] - u_s[i_sol]) / fabs(u_s1[i_sol] + u_s[i_sol]);
		};

		vmax_Iter = std::max_element(errorVector.begin(), errorVector.end());
		SetErrorMax(*vmax_Iter);

		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			u_s[i_sol] = u_s1[i_sol];
		};

		reachedIteration++;
	};

	if (GetIterationsNumber() == GetMaxIterationsNumber())
	{
//		LogConvergenceParameters();
	};

	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = u_s[i_sol];
	};
};

void TBaseGearSolver::UpdateStartingPoints(const TArray u)
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		switch (GetOrder())
		{
		case 2:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u[i_sol]; //new point
			break;
		case 4:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u2[i_sol];
			u2[i_sol] = u3[i_sol];
			u3[i_sol] = u[i_sol]; //new point
			break;
		case 6:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u2[i_sol];
			u2[i_sol] = u3[i_sol];
			u3[i_sol] = u2[i_sol];
			u4[i_sol] = u3[i_sol];
			u5[i_sol] = u[i_sol]; //new point
			break;
		default:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u2[i_sol];
			u2[i_sol] = u3[i_sol];
			u3[i_sol] = u[i_sol]; //new point
			break;
		};
	};
};