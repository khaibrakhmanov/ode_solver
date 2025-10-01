#include "test4.h"
#include "ctest.h"
#include <cmath>
namespace test4
{
	double dydt(double t, std::vector<double> y);
	double exact_sol(double t);
	auto f0 = 0.0;
	auto a = 6.031;
	auto b = 19.74;
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
	auto tau = 1 / b;

	// solver params

	tol = 1e-4;
	dt = 1e-4 * tau;
	dt_min = 1e-6 * tau;
	dt_max = 1e-1 * tau;

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