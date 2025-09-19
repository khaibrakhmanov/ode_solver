#include "test1.h"
#include "ctest.h"

namespace test1
{
	double dydt(double t, std::vector<double> y);
	double exact_sol(double t);
	double f0 = 1.0;
	double alpha = 1.0;
}

CTest1::CTest1()
{
	name = "test1";
}

CTest1::~CTest1()
{
}

void CTest1::SetEqs()
{
	using namespace test1;

	eqs_number = 1;

	PFunctionVectArg rhs = dydt;
	PFunction pSol = &exact_sol;
	ode.AddRHS(rhs);
	ode.AddIC(f0);
	ode.AddSolution(pSol);
}

void CTest1::SetParams()
{
	using namespace test1;

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

double test1::dydt(double t, std::vector<double> y)
{
	return -alpha * y[0];
};

double test1::exact_sol(double t)
{
	return exp(-alpha * t);
};