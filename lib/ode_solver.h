#ifndef ODE_INCLUDED
#define ODE_INCLUDED

#include <vector>

/// @brief user type: one-dimensional vector of real values
typedef std::vector<double> TArray;

/// @brief @brief user type: two-dimensional vector of real values
typedef std::vector<TArray> TArray2D;

/// @brief class: switcher between states 'on/off'
class TSwitcher
{
private:
	/// @brief switcher value
	bool value;
public:
	TSwitcher();

	/// @brief turn on the switcher
	void TurnOn();
	/// @brief tirn off the switcher
	void TurnOff();

	/// @brief Is the switcher on?
	/// @return switcher status (true: on, false: off)
	bool Active();
};

class TBaseODESolver
{
private:

	TSwitcher IncreaseNextStep;
	TSwitcher StepSizeControl;

	double min_step;
	double max_step;

	bool solutionFound;

protected:
	// order of aproximation
	std::size_t order;
	TSwitcher richardsonExtrapolation;

	// solution vestors
	TArray sol; //obtained using one step
	TArray sol2; //obtained using two steps, h/2 and h/2
	TArray solHalf; //at the step h / 2
	// initial condition
	TArray initialVector;
	// errors vector
	TArray errorVector;
	TArray absoluteErrorVector;

	std::size_t variablesNumber;
	// current number of iterations
	std::size_t reachedIteration;
	std::size_t maxIterationsNumber;
	// maximal error among the errors for each component of the solution vector
	double maximalError;

	double dt;
	double time;
	// required tolerance
	double tol;


//	//------------- log -----------------
//	// Log-file ID
//	FILE* pOdeLogFile;
//	/* returns ID of the log file */
//	FILE* Get_file_ID();
//	//-----------------------------------


	void AdjustSolution();

public:
	TBaseODESolver();

	void Complete();

	virtual std::size_t GetOrder() = 0;
	virtual void SetOrder(const std::size_t) = 0;

	bool Solve();

	virtual void AllocateCoefficients() = 0;
	virtual void Step() = 0; //with h
	virtual void Step2() = 0;// with h/2 + h/2

	void ActivateRichardsonExtrapolation();

	void AllocateSolution();

	void Configure();
	void ConfigureDefault();

	void ActivateStepSizeControl();
	bool IfStepSizeControl();

	double GetTimeStep();
	double GetTime();
	void SetTime(const double _time);
	void SetTimeStep(const double _dt);

	void SetTolerance(const double);
	double GetTolerance();
	std::size_t GetIterationsNumber();
	std::size_t GetMaxIterationsNumber();

	/* increase step by factor 2 */
	void DoubleStep();
	/* reduce step by factor 2 */
	void HalveStep();
	void UpdateStep(const double);

	void ReadInitialVector(const TArray, const std::size_t vsize);
	//void ReadInitialVector(const double);
	void CheckSize(const std::size_t _size);
	double GetSolution(const std::size_t);


	/* returns error for solution with given index */
	double GetError(const std::size_t);
	void ComputeError();
	void SetErrorMax(const double err_max);
	double GetErrorMax();

	virtual void UserCheck() = 0;

	std::size_t GetVariablesNumber();
	void SetVariablesNumber(const std::size_t);

	double GetMinStep();
	void SetMinStep(const double);

	double GetMaxStep();
	void SetMaxStep(const double);

	/* returns true if the solution is found, false otherwise */
	bool GetSolutionStatus();


//	//--------- log ------------------
//
//	/* create text file for the errors, messages and warnings reports */
//	void CreateLog(char *fname);
//	/* write title of the current section to the log file */
//	void UpdateLogSection(char *string);
//	/* close log file */
//	void CloseLog();
//	void LogConvergenceParameters();
//	//--------------------------------

};

class TBaseEulerSolver : public TBaseODESolver
{
private:

public:
	TBaseEulerSolver();

	virtual double rhs(const std::size_t _isol, const double t, const TArray y) = 0;

	void Step(); //with h
	void Step2(){};// with h/2 + h/2

	void AllocateCoefficients(){};
	std::size_t GetOrder(){ return 0; };
	void SetOrder(const std::size_t){};
};

class TBasicRungeKuttaSolver : public TBaseODESolver
{
protected:
	//size_t order;
	TArray2D coeff;
	/// y1 = y0 + 0.5 * C0
	TArray y1;
	/// y2 = y1 + 0.5 * C1
	TArray y2;
	/// y3 = y2 + 0.5 * C2
	TArray y3;

	void CalcIntemediateVariables(const TArray y0);

public:
	TBasicRungeKuttaSolver();

	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

	virtual void Step() = 0; //with h
	virtual void Step2() = 0;// with h/2 + h/2

	double GetC(const std::size_t i, const std::size_t j);
	void SetC(const std::size_t i, const std::size_t j, const double arg);

	/// <summary>
	/// returns pointer to the array of a scheme coefficients
	/// </summary>
	/// <returns>pointer to coeff</returns>
	TArray2D* GetCoefficientsPtr();

	void AllocateCoefficients();

	std::size_t GetOrder();
	void SetOrder(const std::size_t);
};

class TBaseRungeKuttaSolver : public TBasicRungeKuttaSolver
{
public:
	TBaseRungeKuttaSolver();

	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

	void Step(); //with h
	void Step2();// with h/2 + h/2
};

class TBaseRungeKuttaFehlberg45Solver : public TBasicRungeKuttaSolver
{
protected:
	// for 4-order step
	double b1_4, b3_4, b4_4, b5_4;
	// for 5-order step
	double b1_5, b3_5, b4_5, b5_5, b6_5;

public:
	TBaseRungeKuttaFehlberg45Solver();

	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

	void Step(); //with h
	void Step2();// with h/2 + h/2
};

class TBaseGearSolver : public TBaseODESolver
{
private:
	//size_t order;
	TArray u0, u1, u2, u3, u4, u5;
	TArray2D u_start;
	TArray u_s, u_s1, u_guess;
	double c_u0, c_u1, c_u2, c_u3, c_u4, c_u5, c_f;
	//TArray error;

	//size_t_vector Iteration;
	//size_t Iteration;

public:
	TBaseGearSolver();

	virtual double rhs(const std::size_t _isol, const double t, const TArray u) = 0;

	void AllocateCoefficients();
	void ReadStartingPoints(const std::size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3, const double _u4, const double _u5);
	void ReadStartingPoints(const std::size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3);
	void ReadStartingPoints(const std::size_t i_sol, const double _u0, const double _u1);
	void ReadStartingPoints(const TArray2D u_in);
	void ReadGuess(const std::size_t i_sol, const double u_g);
	void UpdateStartingPoints(const TArray u);

	double GetExplicitPart(const std::size_t i_sol);
	double GetPhi(const std::size_t i_sol, const TArray u);

	void Step();

	void Step2(){ ; };
	void SetOrder(const std::size_t);
	std::size_t GetOrder();
	//size_t GetIterationsNumber(const size_t i_sol);
};

#endif
