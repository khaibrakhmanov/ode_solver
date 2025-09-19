#include "test5.h"
#include "ctest.h"

namespace test5
{
	double dydt(double t, std::vector<double> y);
	double exact_sol(double t);
	auto f0 = 0.0;
	auto a = 0.026;
	auto b = 4.695e-12;
	auto tau = 0.5 / sqrt(a * b);
}

CTest5::CTest5()
{
	name = "test5";
}

CTest5::~CTest5()
{
}

void CTest5::SetEqs()
{
	using namespace test5;

	eqs_number = 1;

	PFunctionVectArg rhs = dydt;
	PFunction pSol = &exact_sol;
	ode.AddRHS(rhs);
	ode.AddIC(f0);
	ode.AddSolution(pSol);
}

void CTest5::SetParams()
{
	using namespace test5;

	// equation params

	a = 0.026;
	b = 4.695e-12;
	tau = 0.5 / sqrt(a * b);

	// test params
	ic.push_back(f0);
	t_stop = 5e7;

	// solver params

	tol = 1e-4;
	dt = 1e-4 * tau;
	dt_min = 1e-6 * tau;
	dt_max = 1e-1 * tau;

}

double test5::dydt(double t, std::vector<double> y)
{
	return a - b * y[0] * y[0];
};

double test5::exact_sol(double t)
{
	return sqrt(a / b) * (1 - exp(-t / tau)) / (1 + exp(-t / tau));
};