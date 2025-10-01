#pragma once
#include "../lib/rk_solver.h"
#include <fstream>
#include <string>

/// @brief Class - code for testing numerical solution of an ODEs system
class CTest
{
public:
	CTest();
	~CTest();

	/// @brief Run test
	void Run();

protected:
	/// @brief test name
	std::string name;

	/// @brief System of ODEs to be solved
	TODEs ode;
	/// @brief Runge-Kutta solver to solve the system
	TRungeKuttaSolver solver;
	/// @brief Vector of the solution at a given time step
	std::vector<double> sol;
	/// @brief Vector of initial conditions
	std::vector<double> ic;
	/// @brief current time
	double t;
	/// @brief current time step
	double dt;
	/// @brief minimum allowed time step
	double dt_min;
	/// @brief maximum allowed time step
	double dt_max;
	/// @brief time to stop the simulation at
	double t_stop;
	/// @brief required relative accuracy of the integration
	double tol;
	/// @brief number of the equations in the system
	size_t eqs_number;
	/// @brief time step count
	size_t i_t;
	/// @brief name of the file to save the simulation results
	std::ofstream out_file;
	/// @brief variable to store the status of the solution (found or not)
	bool solutionFound;
	/// @brief number of time steps after which the solution is stored on the disc
	size_t savePeriod;

	/// @brief configure equation for a specific test
	virtual void SetEqs() = 0;
	/// @brief configure parameters for a specific test
	virtual void SetParams() = 0;
	/// @brief allocate memory
	void Allocate();
	/// @brief prepare algorythms
	void Prepare();
	/// @brief set the initial conditions
	void SetIC();
	/// @brief execute the solution algorythm
	void Solve();
	/// @brief open a text file to store the solution in
	void OpenResultsFile();
	/// @brief save simulation results at a given time step
	void SaveResults();
};
