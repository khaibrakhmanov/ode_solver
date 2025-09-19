#include "ctest.h"

CTest::CTest()
{

}

CTest::~CTest()
{
}

void CTest::Allocate()
{
	solver.AllocateSolution();
	solver.AllocateCoefficients();
	sol.resize(eqs_number);
}

void CTest::Prepare()
{
	solver.ActivateStepSizeControl();
	solver.SetVariablesNumber(eqs_number);
	solver.SetOrder(4);
	solver.SetTolerance(tol);
	solver.SetMaxStep(dt_max);
	solver.SetMinStep(dt_min);
	solver.ImportODEs(&ode);
}

void CTest::SetIC()
{
	this->t = 0.0;
	i_t = 0;
	solver.SetTimeStep(dt);

	for (size_t i = 0; i < eqs_number; i++)
		sol.at(i) = ic.at(i);
}

void CTest::Solve()
{
	solutionFound = true;
	savePeriod = 1; // to set params!
	OpenResultsFile(); // to prepare!

	do
	{
		solver.ReadInitialVector(sol, eqs_number);
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

		for (size_t i_sol = 0; i_sol < eqs_number; i_sol++)
			sol.at(i_sol) = solver.GetSolution(i_sol);

		if (i_t % 1 == 100)
		{
			printf("time = %5.2e [t0], dt = %5.2e [t0]\n", t, dt);
			printf("  x  = %5.2e [v0]\n", sol[0]);
		};

		/// save current state
		if ((i_t % savePeriod) == 0)
			SaveResults();

		i_t++;
	} while (t <= t_stop);
}

void CTest::OpenResultsFile()
{
	out_file.open(name + "_data.dat", std::ofstream::out);
	if (!out_file.is_open()) {
		throw std::runtime_error("Could not open output data file!");
	}
}

void CTest::SaveResults()
{
	out_file << t << " ";
	for (size_t i_sol = 0; i_sol < eqs_number; i_sol++)
		out_file << sol.at(i_sol) << " ";

	out_file << std::endl;
}

void CTest::Run()
{
	SetEqs();
	SetParams();
	Prepare();
	Allocate();
	SetIC();
	Solve();
}
