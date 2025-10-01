#include <iostream>
#include <algorithm>
#include <cmath>

#include "ode_solver.h"
#include "utypes.h"

TSwitcher::TSwitcher()
{
	value = false;
};

void TSwitcher::TurnOn()
{
	value = true;
};

bool TSwitcher::Active()
{
	return value;
}

void TSwitcher::TurnOff()
{
	value = false;
};

//------------------------------------------------

//FILE* TBaseODESolver::Get_file_ID()
//{
//	return pOdeLogFile;
//}

void TBaseODESolver::AdjustSolution()
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = sol2[i_sol] + absoluteErrorVector[i_sol];
	};
}

TBaseODESolver::TBaseODESolver()
{
	min_step = 1.0e-10;
	max_step = 0.1;
	solutionFound = true;
	maxIterationsNumber = 20;
}

void TBaseODESolver::Complete()
{
//	CloseLog();
};

void TBaseODESolver::ActivateStepSizeControl()
{
	StepSizeControl.TurnOn();
};

bool TBaseODESolver::IfStepSizeControl()
{
	return StepSizeControl.Active();
};

double TBaseODESolver::GetMinStep()
{
	return min_step;
};

double TBaseODESolver::GetMaxStep()
{
	return max_step;
};

void TBaseODESolver::SetMinStep(const double arg)
{
	min_step = arg;
};

void TBaseODESolver::SetMaxStep(const double arg)
{
	max_step = arg;
}

bool TBaseODESolver::GetSolutionStatus()
{
	return solutionFound;
}

//void TBaseODESolver::CreateLog(char * fname)
//{
//	//fopen_s(&pOdeLogFile, fname, "w");
//	pOdeLogFile = fopen(fname, "w");
//
//	if (!pOdeLogFile)
//	{
//		std::cout << "unable to open log file for writing!" << std::endl;
//	};
//}
//
//void TBaseODESolver::UpdateLogSection(char * string)
//{
//	fprintf(Get_file_ID(), "%s\n\n", string);
//}
//
//void TBaseODESolver::CloseLog()
//{
//	fclose(Get_file_ID());
//}
//
//void TBaseODESolver::LogConvergenceParameters()
//{
//	fprintf(Get_file_ID(), "Convergence parameters at t = %5.2e):\n", GetTime());
//	fprintf(Get_file_ID(), "  dt     = %5.2e):\n", GetTimeStep());
//	fprintf(Get_file_ID(), "  eps    = %5.2e):\n", GetErrorMax());
//	fprintf(Get_file_ID(), "  N_iter = %5.2e):\n", GetIterationsNumber());
//
//	fprintf(Get_file_ID(), "-------------------------------------\n");
//};

void TBaseODESolver::SetTime(const double _time)
{
	time = _time;
};

double TBaseODESolver::GetTimeStep()
{
	return dt;
};

void TBaseODESolver::DoubleStep()
{
	SetTimeStep(2.0 * GetTimeStep());
};

void TBaseODESolver::HalveStep()
{
	SetTimeStep(0.5 * GetTimeStep());
}

double TBaseODESolver::GetTime()
{
	return time;
};

void TBaseODESolver::SetTimeStep(const double _dt)
{
	dt = _dt;
};

void TBaseODESolver::SetTolerance(const double _tol)
{
	tol = _tol;
};

void TBaseODESolver::ComputeError()
{
	/* iterator pointing to the maximal value in the vector */
	TArray::iterator vector_max_Iter;
#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(errorVector, sol, sol2)
	{
#pragma omp for nowait
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			errorVector[i_sol] = fabs((sol[i_sol] - sol2[i_sol]) / (fabs(sol2[i_sol]) + fabs(sol[i_sol])));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		//absoluteErrorVector[i_sol] = (sol2[i_sol] - sol[i_sol]) / (pow(2.0, GetOrder()) - 1.0);
		//errorVector[i_sol] = fabs((sol2[i_sol] - sol[i_sol]) / std::max(sol[i_sol], sol2[i_sol]));
		errorVector[i_sol] = fabs((sol2[i_sol] - sol[i_sol]) / sol[i_sol]);
		absoluteErrorVector[i_sol] = (sol2[i_sol] - sol[i_sol]);
	};
#endif

	vector_max_Iter = std::max_element(errorVector.begin(), errorVector.end());
	SetErrorMax(*vector_max_Iter);
}

void TBaseODESolver::SetErrorMax(const double err_max)
{
	maximalError = err_max;
}

double TBaseODESolver::GetErrorMax()
{
	return maximalError;
};

void TBaseODESolver::Configure()
{
	ConfigureDefault();
};

void TBaseODESolver::AllocateSolution()
{
	AllocateVector1D(&initialVector, GetVariablesNumber());
	AllocateVector1D(&sol, GetVariablesNumber());
	AllocateVector1D(&sol2, GetVariablesNumber());
	AllocateVector1D(&solHalf, GetVariablesNumber());
	AllocateVector1D(&errorVector, GetVariablesNumber());
	AllocateVector1D(&absoluteErrorVector, GetVariablesNumber());
};

void TBaseODESolver::CheckSize(const size_t _size)
{
	if (_size != GetVariablesNumber())
		std::cout << "[warning] size of the input vector does not match size of the solution vector" << std::endl;
};

void TBaseODESolver::ReadInitialVector(const TArray _inV, const size_t vsize)
{
	CheckSize(vsize);
	initialVector = _inV;
};

double TBaseODESolver::GetSolution(const size_t index)
{
	return sol[index];
};

double TBaseODESolver::GetError(const size_t index)
{
	return errorVector[index];
};

double TBaseODESolver::GetTolerance()
{
	return tol;
};

bool TBaseODESolver::Solve()
{
	if (!StepSizeControl.Active())
	{
		/* use fixed step size */
		Step();
	}
	else
	{
		/* use automatic step size control */
		Step();
		Step2();
		ComputeError();
		UserCheck();

		reachedIteration = 0;

		if (GetErrorMax() <= GetTolerance())
		{
			if (GetTimeStep() < GetMaxStep())
			{
				IncreaseNextStep.TurnOn();
				if (richardsonExtrapolation.Active())
				{
					AdjustSolution();
				}
			}
		}
		else
		{
			do
			{
				if (GetTimeStep() < GetMinStep())
				{
					solutionFound = false;
					break;
				}
				else
				{
					SetTimeStep(GetTimeStep() * 0.5);
					Step();
					Step2();
					ComputeError();
					reachedIteration++;
				}
			} while ((GetErrorMax() > GetTolerance()) && (GetIterationsNumber() < GetMaxIterationsNumber()));
			if ((GetIterationsNumber() < GetMaxIterationsNumber()) && richardsonExtrapolation.Active())
				AdjustSolution();
		}

		if (GetIterationsNumber() == GetMaxIterationsNumber())
		{
//			LogConvergenceParameters();
		};

	};

	return GetSolutionStatus();
};


size_t TBaseODESolver::GetIterationsNumber()
{
	return reachedIteration;
}

size_t TBaseODESolver::GetMaxIterationsNumber()
{
	return maxIterationsNumber;
}

size_t TBaseODESolver::GetVariablesNumber()
{
	return variablesNumber;
};

void TBaseODESolver::SetVariablesNumber(const size_t arg)
{
	variablesNumber = arg;
};

void TBaseODESolver::ConfigureDefault()
{
	SetVariablesNumber(variablesNumber);
	IncreaseNextStep.TurnOff();
};

void TBaseODESolver::UpdateStep(const double old_step)
{
	if (IncreaseNextStep.Active())
	{
		SetTimeStep(old_step * 2.0);
		IncreaseNextStep.TurnOff();
	};
};

void TBaseODESolver::ActivateRichardsonExtrapolation()
{
	richardsonExtrapolation.TurnOn();
};