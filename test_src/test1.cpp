#include "test1.h"
#include "ctest.h"

double dydt(double t, std::vector<double> y);
double exact_sol(double t);
double f0 = 1.0;
double alpha = 1.0;

void CTest1::SetEqs()
{
	eqs_number = 1;

	PFunctionVectArg rhs = dydt;
	PFunction pSol = &exact_sol;
	ode.AddRHS(rhs);
	ode.AddIC(f0);
	ode.AddSolution(pSol);
}

void CTest1::SetParams()
{
	// equation params

	alpha = 1.0;

	// solver params

	tol = 1e-4;
	dt = 1e-4;
	dt_min = 1e-4;
	dt_max = 1e-2;

	// test params
	ic.push_back(f0);
	t_stop = 100.0;
}

double dydt(double t, std::vector<double> y)
{
	return -alpha * y[0];
};

double exact_sol(double t)
{
	return exp(-alpha * t);
};