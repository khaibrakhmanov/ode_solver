#include "euler_solver.h"

TBaseEulerSolver::TBaseEulerSolver()
{

};

void TBaseEulerSolver::Step()
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() * rhs(i_sol, GetTime(), initialVector);
	};
};