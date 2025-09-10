#include <iostream>
#include <fstream>
#include "solver.h"
#include "eqs.h"

/// @brief Name of the file to save the simulation results
std::ofstream outfile;

const double alpha = 1.0;

/// @brief Open output data file
void OpenOutputDataFile();
/// @brief Save the simulation results at current time step
void SaveCurrentData(const double t, const double sol, const double an_sol);

double f1(double t, std::vector<double> y);
double f2(double t, std::vector<double> y);
double sol1(double t);
double sol2(double t);

int main(int argc, char* argv[])
{
	// if argv are entered then the parameter order is following:
	// argv[1] -> IC, i.e. y(t=0)
	// argv[2] -> alpha, i.e. parameter of the equation dy/dt = -alpha*dt
	// argv[3] -> dt, i.e. (initial) time step

	std::cout << "You have entered " << argc << " arguments:" << std::endl;
	for (size_t i = 0; i < argc; i++)
		std::cout << argv[i] << std::endl;

	std::cout << "Press Enter to start a simulation" << std::endl;
	std::cin.get();

	double y0 = 0.0;
	TExpODE eq1(1000.0, y0);
	TLinODE eq2(6.031, 19.74, y0);
	TQuadODE eq3(0.026, 4.695e-12, y0);
	
	
	TODEs eqs;
	PFunction ode1 = f1;
	PFunction ode2 = f2;
	//PFunction ode1_sol = sol1;

	eqs.AddRHS(ode1);
	eqs.AddIC(1.0);
	eqs.AddSolution(&sol1);

	eqs.AddRHS(ode2);
	eqs.AddIC(1.0);
	eqs.AddSolution(&sol2);

	/// time of the simulation stop
	const auto t_stop = 3e8;

	TRungeKuttaSolver solver;
	std::vector<double> sol;
	size_t sol_size = 1;
	solver.ActivateStepSizeControl();
	solver.SetVariablesNumber(1);
	solver.SetOrder(4);
	solver.SetTolerance(1e-4);
	solver.SetMaxStep(1e-2 * t_stop);
	solver.SetMinStep(1e-4 * t_stop);
	solver.AllocateSolution();
	solver.AllocateCoefficients();
	sol.resize(1);

	solver.ImportODE(&eq3);

	/// time step count
	auto i_t = 1;
	/// current time
	auto t = 0.0;
	sol[0] = y0;
	auto an_sol = y0;
	/// current time step
	auto dt = 1e-4 * t_stop;
	solver.SetTimeStep(dt);

	bool solutionFound = true;
	auto savePeriod = 1;
	OpenOutputDataFile();

	do
	{
		solver.ReadInitialVector(sol, sol_size);
		solutionFound = solver.Solve();
		if (!solutionFound)
		{
			printf("Warning! Solution was not found. Exiting.\n");
			break;
		}
		dt = solver.GetTimeStep();
		solver.UpdateStep(dt);

		t += dt;
		solver.SetTime(t);

		sol[0] = solver.GetSolution(0);
		an_sol = eq3.Solution(t);

		if (i_t % 1 == 100)
		{
			printf("time = %5.2e [t0], dt = %5.2e [t0]\n", t, dt);
			printf("  x  = %5.2e [v0]\n", sol[0]);
		};

		/// save current state
		if ((i_t % savePeriod) == 0)
			SaveCurrentData(t, sol[0], an_sol);

		i_t++;
	} while (t <= t_stop);

	return 0;
}

void OpenOutputDataFile()
{
    outfile.open("data.dat", std::ofstream::out);
    if (!outfile.is_open()) {
        throw std::runtime_error("Could not open output data file!");
    }
}

void SaveCurrentData(const double t, const double sol, const double an_sol)
{
    outfile << t << " " << sol << " " << an_sol << std::endl;
}

double f1(double t, std::vector<double> y)
{
	return -alpha * y.at(1);
};

double f2(double t, std::vector<double> y)
{
	return alpha * y.at(0) - y.at(1);
};

double sol1(double t)
{
	double b = std::sqrt(4.0 * alpha * alpha - 1.0);
	return std::exp(-0.5 * t) * ((1.0 - 2.0 * alpha) * std::sin(0.5 * b * t) / b + std::cos(0.5 * b * t));
};

double sol2(double t)
{
	double b = std::sqrt(4.0 * alpha * alpha - 1.0);
	return std::exp(-0.5 * t) * ((2.0 * alpha - 1.0) * std::sin(0.5 * b * t) / b + std::cos(0.5 * b * t));
};