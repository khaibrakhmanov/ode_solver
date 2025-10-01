#include "test2.h"
#include "ctest.h"
#include <cmath>
namespace test2
{
	double dydt_1(double t, std::vector<double> y);
	double dydt_2(double t, std::vector<double> y);
	double exact_sol_1(double t);
	double exact_sol_2(double t);
	double f0_1 = 1.0;
	double f0_2 = 1.0;
	double alpha = 1000.0;
}

CTest2::CTest2()
{
	name = "test2";
}

CTest2::~CTest2()
{
}

void CTest2::SetEqs()
{
	using namespace test2;

	eqs_number = 2;

	PFunctionVectArg rhs_1 = dydt_1;
	PFunctionVectArg rhs_2 = dydt_2;
	PFunction pExactSol_1 = &exact_sol_1;
	PFunction pExactSol_2 = &exact_sol_2;

	ode.AddRHS(rhs_1);
	ode.AddIC(f0_1);
	ode.AddSolution(pExactSol_1);

	ode.AddRHS(rhs_2);
	ode.AddIC(f0_2);
	ode.AddSolution(pExactSol_2);
}

void CTest2::SetParams()
{
	using namespace test2;

	// equation params

	alpha = 1000.0;
	ode.SetParameter("alpha", alpha);

	// solver params
	auto tau = 1 / alpha;

	tol = 1e-4;
	dt = 1e-4 * tau;
	dt_min = 1e-6 * tau;
	dt_max = 1 * tau;

	// test params
	ic.push_back(f0_1);
	ic.push_back(f0_2);
	t_stop = 0.02;
}

double test2::dydt_1(double t, std::vector<double> y)
{
	return -alpha * y[0];
};

double test2::dydt_2(double t, std::vector<double> y)
{
	return -y[1];
};

double test2::exact_sol_1(double t)
{
	return f0_1 * exp(-alpha * t);
};

double test2::exact_sol_2(double t)
{
	return f0_2 * exp(-t);
};