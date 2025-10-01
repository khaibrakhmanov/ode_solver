#include "test3.h"
#include "ctest.h"
#include <cmath>
namespace test3
{
	double dydt_1(double t, std::vector<double> y);
	double dydt_2(double t, std::vector<double> y);
	double exact_sol_1(double t);
	double exact_sol_2(double t);
	auto f0_1 = 1.0;
	auto f0_2 = 1.0;
	auto alpha = 1000.0;
	auto b = sqrt(4.0 * alpha * alpha - 1.0);
}

CTest3::CTest3()
{
	name = "test3";
}

CTest3::~CTest3()
{
}

void CTest3::SetEqs()
{
	using namespace test3;

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

void CTest3::SetParams()
{
	using namespace test3;

	// equation params

	alpha = 1000.0;
	b = sqrt(4.0 * alpha * alpha - 1.0);
	auto tau = 2 / b;
	ode.SetParameter("alpha", alpha);

	// solver params

	tol = 1e-4;
	dt = 1e-4 * tau;
	dt_min = 1e-6 * tau;
	dt_max = 1e-1 * tau;

	// test params
	ic.push_back(f0_1);
	ic.push_back(f0_2);
	t_stop = 0.05;
}

double test3::dydt_1(double t, std::vector<double> y)
{
	return -alpha * y[1];
};

double test3::dydt_2(double t, std::vector<double> y)
{
	return alpha * y[0] - y[1];
};

double test3::exact_sol_1(double t)
{
	auto b = sqrt(4.0 * alpha * alpha - 1.0);
	return exp(-0.5 * t) * ((1.0 - 2.0 * alpha) * sin(0.5 * b * t) / b + cos(0.5 * b * t));
};

double test3::exact_sol_2(double t)
{
	auto b = sqrt(4.0 * alpha * alpha - 1.0);
	return exp(-0.5 * t) * ((2.0 * alpha - 1.0) * sin(0.5 * b * t) / b + cos(0.5 * b * t));
};