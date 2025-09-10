#include <iostream>
#include <algorithm>
#include <cmath>

#include "ode_solver.h"

/// @brief Allocate memory for a 1-dimensional array
/// @param v array to allocate memory for
/// @param N number of elements
void AllocateVector1D(TArray *v, const size_t N);

/// @brief Allocate memory for a 2-dimensional array
/// @param v array to allocate memory for
/// @param N1 number of elements in the 1-st dimension
/// @param N2 number of elements in the 2-nd dimension
void AllocateVector2D(TArray2D *v, const size_t N1, const size_t N2);

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
		errorVector[i_sol] = fabs((sol2[i_sol] - sol[i_sol]) / std::max(sol[i_sol], sol2[i_sol]));
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

//------------------------------------

TBasicRungeKuttaSolver::TBasicRungeKuttaSolver()
{
	richardsonExtrapolation.TurnOff();
}

double TBasicRungeKuttaSolver::GetC(const size_t i, const size_t j)
{
	return coeff[i][j];
};

void TBasicRungeKuttaSolver::SetC(const size_t i, const size_t j, const double arg)
{
	coeff[i][j] = arg;
};

TArray2D* TBasicRungeKuttaSolver::GetCoefficientsPtr()
{
	return &coeff;
}

void TBasicRungeKuttaSolver::AllocateCoefficients()
{
	AllocateVector2D(&coeff, GetOrder(), GetVariablesNumber());
};

size_t TBasicRungeKuttaSolver::GetOrder()
{
	return order;
};

void TBasicRungeKuttaSolver::SetOrder(const size_t arg)
{
	order = arg;
};

void TBasicRungeKuttaSolver::CalcIntemediateVariables(const TArray y0)
{
    for(std::size_t i = 0; i < GetVariablesNumber(); i++)
    {
        y1[i] = y0[i] + 0.5 * dt * GetC(0, i);
        y2[i] = y0[i] + 0.5 * dt * GetC(1, i);
        y3[i] = y0[i] + dt * GetC(2, i);
    }
}

//------------- Ordinary RK solver----------------------
TBaseRungeKuttaSolver::TBaseRungeKuttaSolver()
{
	SetOrder(4);
}

void TBaseRungeKuttaSolver::Step()
{
	CalcCoefficients(initialVector, GetTimeStep(), GetTime());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(sol, initialVector)
	{
#pragma omp for nowait
		for (int i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			sol[i_sol] = initialVector[i_sol] + GetTimeStep() / 6.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() / 6.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

};

void TBaseRungeKuttaSolver::Step2()
{

	CalcCoefficients(initialVector, 0.5 * GetTimeStep(), GetTime());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(solHalf, initialVector)
	{
#pragma omp for nowait
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			solHalf[i_sol] = initialVector[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		solHalf[i_sol] = initialVector[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

	CalcCoefficients(solHalf, 0.5 * GetTimeStep(), GetTime() + 0.5 * GetTimeStep());

#ifdef PARALLEL_MODE
#pragma omp parallel num_threads(NUM_THREADS) shared(sol2, solHalf)
	{
#pragma omp for nowait
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			sol2[i_sol] = solHalf[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
		};
	};
#else
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol2[i_sol] = solHalf[i_sol] + GetTimeStep() / 12.0 * (GetC(0, i_sol) + 2.0 * GetC(1, i_sol) + 2.0 * GetC(2, i_sol) + GetC(3, i_sol));
	};
#endif

};

//------------------------------

TBaseEulerSolver::TBaseEulerSolver()
{

};

void TBaseEulerSolver::Step()
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() * rhs(i_sol, GetTime(), initialVector);
	};
};

//------------ Gear method -------------

void TBaseGearSolver::SetOrder(const size_t _order)
{
	order = _order;

	switch (GetOrder())
	{
	case 2:
		c_u0 = 1.0 / 3.0;
		c_u1 = 4.0 / 3.0;
		c_f = 2.0 / 3.0;
		break;
	case 4:
		c_u0 = 3.0 / 25.0;
		c_u1 = 16.0 / 25.0;
		c_u2 = 36.0 / 25.0;
		c_u3 = 48.0 / 25.0;
		c_f = 12.0 / 25.0;
		break;
	case 6:
		c_u0 = 10.0 / 147.0;
		c_u1 = 72.0 / 147.0;
		c_u2 = 225.0 / 147.0;
		c_u3 = 400.0 / 147.0;
		c_u4 = 450.0 / 147.0;
		c_u5 = 360.0 / 147.0;
		c_f = 60.0 / 147.0;
		break;
	default:
		c_u0 = 3.0 / 25.0;
		c_u1 = 16.0 / 25.0;
		c_u2 = 36.0 / 25.0;
		c_u3 = 48.0 / 25.0;
		c_f = 12.0 / 25.0;
		break;
	}
};

size_t TBaseGearSolver::GetOrder()
{
	return order;
};

TBaseGearSolver::TBaseGearSolver()
{
	maxIterationsNumber = 20;
};

void TBaseGearSolver::AllocateCoefficients()
{
	AllocateVector1D(&u0, GetVariablesNumber());
	AllocateVector1D(&u1, GetVariablesNumber());
	AllocateVector1D(&u2, GetVariablesNumber());
	AllocateVector1D(&u3, GetVariablesNumber());
	AllocateVector1D(&u4, GetVariablesNumber());
	AllocateVector1D(&u5, GetVariablesNumber());
	AllocateVector1D(&u_s, GetVariablesNumber());
	AllocateVector1D(&u_s1, GetVariablesNumber());
	AllocateVector1D(&u_guess, GetVariablesNumber());
};

double TBaseGearSolver::GetExplicitPart(const size_t i_sol)
{
	switch (GetOrder())
	{
	case 2:
		return c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	case 4:
		return c_u3 * u3[i_sol] - c_u2 * u2[i_sol] + c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	case 6:
		return c_u5 * u5[i_sol] - c_u4 * u4[i_sol] + c_u3 * u3[i_sol] - c_u2 * u2[i_sol] + c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	default:
		return c_u3 * u3[i_sol] - c_u2 * u2[i_sol] + c_u1 * u1[i_sol] - c_u0 * u0[i_sol];
		break;
	}
};

double TBaseGearSolver::GetPhi(const size_t i_sol, const TArray u)
{
	return GetExplicitPart(i_sol) + c_f * GetTimeStep() * rhs(i_sol, GetTime(), u);
};

void TBaseGearSolver::ReadStartingPoints(const size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3)
{
	u0[i_sol] = _u0;
	u1[i_sol] = _u1;
	u2[i_sol] = _u2;
	u3[i_sol] = _u3;
};

void TBaseGearSolver::ReadStartingPoints(const size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3, const double _u4, const double _u5)
{
	u0[i_sol] = _u0;
	u1[i_sol] = _u1;
	u2[i_sol] = _u2;
	u3[i_sol] = _u3;
	u4[i_sol] = _u4;
	u5[i_sol] = _u5;
};

void TBaseGearSolver::ReadStartingPoints(const size_t i_sol, const double _u0, const double _u1)
{
	u0[i_sol] = _u0;
	u1[i_sol] = _u1;
};

void TBaseGearSolver::ReadStartingPoints(const TArray2D u_in)
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		for (size_t i = 0; i < GetOrder(); i++)
		{
			u_start[i_sol][i] = u_in[i_sol][i];
		};
	};
};

void TBaseGearSolver::ReadGuess(const size_t i_sol, const double u_g)
{
	u_guess[i_sol] = u_g;
}

void TBaseGearSolver::Step()
{
	TArray::iterator vmax_Iter;
	SetErrorMax(1.0);

	//initial iteration
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		switch (GetOrder())
		{
		case 2:
			u_s[i_sol] = u1[i_sol];
			break;
		case 4:
			if (u2[i_sol] > 1e-20)
				u_s[i_sol] = u3[i_sol] * u3[i_sol] / u2[i_sol];
			else
				u_s[i_sol] = 2.0 * u3[i_sol] - u2[i_sol];
			break;
		case 6:
			u_s[i_sol] = u5[i_sol] * u5[i_sol] / u4[i_sol];
			break;
		default:
			u_s[i_sol] = u3[i_sol] * u3[i_sol] / u2[i_sol];
			break;
		};
	};

	reachedIteration = 0;
	//iteration cycle
	while ((GetErrorMax() >= GetTolerance()) && (reachedIteration < maxIterationsNumber))
	{
		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			u_s1[i_sol] = GetPhi(i_sol, u_s);
			errorVector[i_sol] = fabs(u_s1[i_sol] - u_s[i_sol]) / fabs(u_s1[i_sol] + u_s[i_sol]);
		};

		vmax_Iter = std::max_element(errorVector.begin(), errorVector.end());
		SetErrorMax(*vmax_Iter);

		for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
		{
			u_s[i_sol] = u_s1[i_sol];
		};

		reachedIteration++;
	};

	if (GetIterationsNumber() == GetMaxIterationsNumber())
	{
//		LogConvergenceParameters();
	};

	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = u_s[i_sol];
	};
};

void TBaseGearSolver::UpdateStartingPoints(const TArray u)
{
	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		switch (GetOrder())
		{
		case 2:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u[i_sol]; //new point
			break;
		case 4:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u2[i_sol];
			u2[i_sol] = u3[i_sol];
			u3[i_sol] = u[i_sol]; //new point
			break;
		case 6:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u2[i_sol];
			u2[i_sol] = u3[i_sol];
			u3[i_sol] = u2[i_sol];
			u4[i_sol] = u3[i_sol];
			u5[i_sol] = u[i_sol]; //new point
			break;
		default:
			u0[i_sol] = u1[i_sol];
			u1[i_sol] = u2[i_sol];
			u2[i_sol] = u3[i_sol];
			u3[i_sol] = u[i_sol]; //new point
			break;
		};
	};
};

//-------- Runge-Kutta-Fehlberg 4(5) solver -----------------------
TBaseRungeKuttaFehlberg45Solver::TBaseRungeKuttaFehlberg45Solver()
{
	SetOrder(5);
	b1_4 = 25.0 / 216.0;
	b3_4 = 1408.0 / 2565.0;
	b4_4 = 2197.0 / 4164.0;
	b5_4 = -0.2;

	b1_5 = 16.0 / 135.0;
	b3_5 = 6656.0 / 12825.0;
	b4_5 = 28561.0 / 56430.0;
	b5_5 = -9.0 / 50.0;
	b6_5 = 2.0 / 55.0;
}

void TBaseRungeKuttaFehlberg45Solver::Step()
{
	CalcCoefficients(initialVector, GetTimeStep(), GetTime());

	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol[i_sol] = initialVector[i_sol] + GetTimeStep() * (b1_4 * GetC(0, i_sol) + b3_4 * GetC(2, i_sol) + b4_4 * GetC(3, i_sol) - b5_4 * GetC(4, i_sol));
	};
}

void TBaseRungeKuttaFehlberg45Solver::Step2()
{
	CalcCoefficients(initialVector, GetTimeStep(), GetTime());

	for (size_t i_sol = 0; i_sol < GetVariablesNumber(); i_sol++)
	{
		sol2[i_sol] = initialVector[i_sol] + GetTimeStep() * (b1_5 * GetC(0, i_sol) + b3_5 * GetC(2, i_sol) + b4_5 * GetC(3, i_sol) + b5_5 * GetC(4, i_sol) + b6_5 * GetC(5, i_sol));
	};
}

// -----------------------------------------------------------------------------------------

void AllocateVector1D(TArray *v, const size_t N)
{
	(*v).resize(N);
}

void AllocateVector2D(TArray2D *v, const size_t N1, const size_t N2) {
	(*v).resize(N1);
	for (size_t i = 0; i < N1; i++)
	{
		(*v)[i].resize(N2);
	};
};
