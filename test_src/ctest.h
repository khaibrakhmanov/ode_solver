#pragma once
#include "eqs.h"
#include "solver.h"
#include <fstream>

class CTest
{
public:
	CTest();
	~CTest();

	void Run();

protected:
	TODEs ode;
	TRungeKuttaSolverSystem solver;
	std::vector<double> sol;
	std::vector<double> ic;
	/// current time
	double t, dt, dt_min, dt_max, t_stop;
	double tol;
	size_t eqs_number;
	/// time step count
	size_t i_t;
	/// @brief Name of the file to save the simulation results
	std::ofstream out_file;
	bool solutionFound;
	size_t savePeriod;

	virtual void SetEqs() = 0;
	virtual void SetParams() = 0;
	void Allocate();
	void Prepare();
	void SetIC();
	void Solve();
	void OpenResultsFile();
	void SaveResults();

};
