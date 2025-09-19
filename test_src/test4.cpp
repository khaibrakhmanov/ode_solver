#include "test4.h"
#include "ctest.h"

namespace test4
{
	double dydt(double t, std::vector<double> y);
	double exact_sol(double t);
	double f0 = 0.0;
	double a = 6.031;
	double b = 19.74;
}

CTest4::CTest4()
{
	name = "test4";
}

CTest4::~CTest4()
{
}

void CTest4::SetEqs()
{
	using namespace test4;

	eqs_number = 1;

	PFunctionVectArg rhs = dydt;
	PFunction pSol = &exact_sol;
	ode.AddRHS(rhs);
	ode.AddIC(f0);
	ode.AddSolution(pSol);
}

void CTest4::SetParams()
{
	using namespace test4;

	// equation params

	a = 6.031;
	b = 19.74;

	// solver params

	tol = 1e-4;
	dt = 1e-4;
	dt_min = 1e-4;
	dt_max = 1e-2;

	// test params
	ic.push_back(f0);
	t_stop = 1.0;
}

double test4::dydt(double t, std::vector<double> y)
{
	return a - b * y[0];
};

double test4::exact_sol(double t)
{
	const auto tau = 1 / b;
	return (a / b) * (1 - exp(-t / tau));
};