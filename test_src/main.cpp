#include <iostream>
#include <fstream>
#include "solver.h"
#include "eqs.h"

/// @brief Name of the file to save the simulation results
std::ofstream outfile;

/// @brief Open output data file
void OpenOutputDataFile();
/// @brief Save the simulation results at current time step
void SaveCurrentData(const double t, const double sol);


int main()
{
	double x0 = 1.0;
	TExpODE eq1(1.0, 1.0);
	TRungeKuttaSolver solver;
	std::vector<double> sol;
	size_t sol_size = 1;

	//solver.ActivateStepSizeControl();
	solver.SetVariablesNumber(1);
	solver.SetOrder(4);
	solver.SetTolerance(1e-4);
	solver.SetMaxStep(0.1);
	solver.SetMinStep(1e-4);
	solver.AllocateSolution();
	solver.AllocateCoefficients();
	sol.resize(1);

	solver.ImportODE(&eq1);

	/// time step count
	auto i_t = 1;
	/// current time
	auto t = 0.0;
	sol[0] = x0;
	/// current time step
	auto dt = 1e-4;
	solver.SetTimeStep(dt);
	/// time of the simulation stop
	const auto t_stop = 1.0;

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

		if (i_t % 1 == 0)
		{
			printf("time = %5.2e [t0], dt = %5.2e [t0]\n", t, dt);
			printf("  x  = %5.2e [v0]\n", sol[0]);
		};

		/// save current state
		if ((i_t % savePeriod) == 0)
			SaveCurrentData(t, sol[0]);

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

void SaveCurrentData(const double t, const double sol)
{
    outfile << t << " " << sol << std::endl;
}
