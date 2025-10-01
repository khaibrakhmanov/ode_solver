#ifndef ODE_INCLUDED
#define ODE_INCLUDED

#include "utypes.h"

/// @brief class: switcher between states 'on/off'
class TSwitcher
{
private:
	/// @brief switcher value
	bool value;
public:
	TSwitcher();

	/// @brief Turn on the switcher
	void TurnOn();
	/// @brief Torn off the switcher
	void TurnOff();

	/// @brief Is the switcher on?
	/// @return switcher status (true: on, false: off)
	bool Active();
};

/// @brief Class: basic ODE solver
class TBaseODESolver
{
private:
	/// @brief switcher: increase next time step (on) or not (off)?
	TSwitcher IncreaseNextStep;
	/// @brief switcher: use step size control (on) or not (off)?
	TSwitcher StepSizeControl;
	/// @brief minimum allowed time step
	double min_step;
	/// @brief maximum allowed time step
	double max_step;
	/// @brief solution status: found (true) or not (false)?
	bool solutionFound;

protected:
	// order of the aproximation of the scheme
	std::size_t order;
	/// @brief switcher: use Richardson extrapolation (on) or not (off)?
	TSwitcher richardsonExtrapolation;

	// --- solution vectors -----

	TArray sol;     // solution obtained using one step dt = h
	TArray sol2;    // solution obtained using two steps, dt1=h/2 and dt2=h/2
	TArray solHalf; // solution at the step dt = h / 2

	// --------------------------
	
	// vector of the initial conditions
	TArray initialVector;
	// vector of relative errors
	TArray errorVector;
	// vector of absolute errors
	TArray absoluteErrorVector;

	/// @brief number of variables in the system of equations to be solved
	std::size_t variablesNumber;
	// current number of iterations
	std::size_t reachedIteration;
	/// @brief maximum allowed number of iterations
	std::size_t maxIterationsNumber;
	/// @brief maximal error among the errors for each component of the solution vector
	double maximalError;
	/// @brief current time step
	double dt;
	/// @brief current time
	double time;
	/// @brief required tolerance (accuracy of the calculation) in the scheme with step size control
	double tol;


//	//------------- log -----------------
//	// Log-file ID
//	FILE* pOdeLogFile;
//	/* returns ID of the log file */
//	FILE* Get_file_ID();
//	//-----------------------------------

	/// @brief Adjuct the solution if the Richardson extrapolation is used
	void AdjustSolution();

public:
	TBaseODESolver();

	/// @brief Finalise the solver execution
	void Complete();

	/// @brief Returns the order of the scheme used
	/// @return order
	virtual std::size_t GetOrder() = 0;

	/// @brief Set the order of the scheme to be used
	/// @param order
	virtual void SetOrder(const std::size_t) = 0;

	/// @brief Main algorythm: solution of the equations using corresponding scheme
	/// @return Solution status (true: solution is found, false: not)
	bool Solve();

	/// @brief Allocate internal memory for coefficients of the scheme
	virtual void AllocateCoefficients() = 0;

	/// @brief Make one step of integration with the current value of time step (dt)
	virtual void Step() = 0;

	/// @brief Make one step of integration using two substeps (dt = 0.5dt + 0.5dt)
	virtual void Step2() = 0;

	/// @brief Activate Richardson extrapolation
	void ActivateRichardsonExtrapolation();

	/// @brief Allocate memory for internal solution vectors
	void AllocateSolution();

	/// @brief Configure solver
	void Configure();

	/// @brief Set default solver settings when configuring
	void ConfigureDefault();

	/// @brief Activate automatic control of the time step size
	void ActivateStepSizeControl();

	/// @brief Check, if step size control is activated
	/// @return True: step size control is On, False: Off
	bool IfStepSizeControl();

	/// @brief Returns current time step
	/// @return dt
	double GetTimeStep();

	/// @brief Returns current time
	/// @return t
	double GetTime();

	/// @brief Set current time from input
	/// @param _time time to be set
	void SetTime(const double _time);

	/// @brief Set time step from input
	/// @param _dt time step
	void SetTimeStep(const double _dt);

	/// @brief Set required relative accuracy from input
	/// @param  relative accuracy (tolerance)
	void SetTolerance(const double);

	/// @brief Returns relative accuracy used by the solver
	/// @return tol
	double GetTolerance();

	/// @brief Returns number of current (time step subcycling) iteration reached
	/// @return reachedIteration
	std::size_t GetIterationsNumber();

	/// @brief Returns maximum allowed number of iterations
	/// @return maxIterationsNumber
	std::size_t GetMaxIterationsNumber();

	/// @brief Increase time step by factor 2
	void DoubleStep();

	/// @brief Reduce time step by factor 2
	void HalveStep();

	/// @brief Update current value of time step from input
	/// @param dt time step to be used
	void UpdateStep(const double dt);

	/// @brief Read vector of initial values from input
	/// @param v0 vector of initial values
	/// @param vsize vector size
	void ReadInitialVector(const TArray v0, const std::size_t vsize);

	/// @brief Check if the size agrees with the number of variables
	/// @param _size size to be checked
	void CheckSize(const std::size_t _size);

	/// @brief Returns the solution of a given equation index
	/// @param i_sol equation index
	/// @return solution[i_sol]
	double GetSolution(const std::size_t i_sol);

	/// @brief Returns error for the solution with given index (in a calculation with step size control)
	/// @param i_sol equation index
	/// @return error[i_sol]
	double GetError(const std::size_t i_sol);

	/// @brief Compute the error vectors at the current time step
	void ComputeError();

	/// @brief Set maximum allowed value of the error in time step subcycling
	/// @param err_max maximum error
	void SetErrorMax(const double err_max);

	/// @brief Get maximum value of the error in time step subcycling
	/// @return maximal error allowed
	double GetErrorMax();

	/// @brief Some user defined procedures before the start of the simulation
	virtual void UserCheck() = 0;

	/// @brief Get the number of variables in the system to be solved
	/// @return N
	std::size_t GetVariablesNumber();

	/// @brief Set the number of variables in the system to be solved from the input
	/// @param n_var number of variables
	void SetVariablesNumber(const std::size_t n_var);

	/// @brief Get minimum allowed time step in the solver with step size control
	/// @return dt_min
	double GetMinStep();

	/// @brief Set minimum allowed time step in the solver with step size control
	/// @param min_step minimum time step
	void SetMinStep(const double min_step);

	/// @brief Get maximum allowed time step in the solver with step size control
	/// @return dt_max
	double GetMaxStep();

	/// @brief Set maximum allowed time step in the solver with step size control
	/// @param max_step maximum time step
	void SetMaxStep(const double max_step);

	/// @brief Check the solution status
	/// @return True: solution is found, False: not
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

// /// @brief Class: basic Euler solver
// class TBaseEulerSolver : public TBaseODESolver
// {
// private:

// public:
// 	TBaseEulerSolver();

// 	virtual double rhs(const std::size_t _isol, const double t, const TArray y) = 0;

// 	void Step(); //with h
// 	void Step2(){};// with h/2 + h/2

// 	void AllocateCoefficients(){};
// 	std::size_t GetOrder(){ return 0; };
// 	void SetOrder(const std::size_t){};
// };

// class TBasicRungeKuttaSolver : public TBaseODESolver
// {
// protected:
// 	/// @brief Matrix of the RK coefficients
// 	TArray2D coeff;

// 	/// intermediate solution: y1 = y0 + 0.5 * C0
// 	TArray y1;
// 	/// intermediate solution: y2 = y1 + 0.5 * C1
// 	TArray y2;
// 	/// intermediate solution: y3 = y2 + 0.5 * C2
// 	TArray y3;
// 	/// @brief Calculate intermediate vectors at current time step
// 	/// @param y0 vector of initial values
// 	void CalcIntemediateVariables(const TArray y0);

// public:
// 	TBasicRungeKuttaSolver();

// 	/// @brief Calculate the coefficients of the scheme
// 	/// @param y0 vector of the initial values
// 	/// @param _dt current time step
// 	/// @param _time surrent time
// 	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

// 	//virtual void Step() = 0; //with h
// 	//virtual void Step2() = 0;// with h/2 + h/2

// 	/// @brief Returns the RK coefficient for given indeces
// 	/// @param i index of the coefficient
// 	/// @param j index of the equation
// 	/// @return coeff[i, j]
// 	double GetC(const std::size_t i, const std::size_t j);

// 	/// @brief Set the RK coefficient for given indeces
// 	/// @param i index of the coefficient
// 	/// @param j index of the equation
// 	/// @param arg value of the coefficient
// 	void SetC(const std::size_t i, const std::size_t j, const double arg);

// 	/// <summary>
// 	/// Returns pointer to the array of the scheme coefficients
// 	/// </summary>
// 	/// <returns>pointer to coeff</returns>
// 	TArray2D* GetCoefficientsPtr();

// 	/// @brief Allocate memory for the scheme coefficients
// 	void AllocateCoefficients();

// 	/// @brief Get the order of the scheme
// 	/// @return order
// 	std::size_t GetOrder();

// 	/// @brief Get the order of the scheme
// 	/// @param o order
// 	void SetOrder(const std::size_t o);
// };

// class TBaseRungeKuttaSolver : public TBasicRungeKuttaSolver
// {
// public:
// 	TBaseRungeKuttaSolver();

// 	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

// 	void Step(); //with h
// 	void Step2();// with h/2 + h/2
// };

// class TBaseRungeKuttaFehlberg45Solver : public TBasicRungeKuttaSolver
// {
// protected:
// 	// for 4-order step
// 	double b1_4, b3_4, b4_4, b5_4;
// 	// for 5-order step
// 	double b1_5, b3_5, b4_5, b5_5, b6_5;

// public:
// 	TBaseRungeKuttaFehlberg45Solver();

// 	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

// 	void Step(); //with h
// 	void Step2();// with h/2 + h/2
// };

// class TBaseGearSolver : public TBaseODESolver
// {
// private:
// 	//size_t order;
// 	TArray u0, u1, u2, u3, u4, u5;
// 	TArray2D u_start;
// 	TArray u_s, u_s1, u_guess;
// 	double c_u0, c_u1, c_u2, c_u3, c_u4, c_u5, c_f;
// 	//TArray error;

// 	//size_t_vector Iteration;
// 	//size_t Iteration;

// public:
// 	TBaseGearSolver();

// 	virtual double rhs(const std::size_t _isol, const double t, const TArray u) = 0;

// 	void AllocateCoefficients();
// 	void ReadStartingPoints(const std::size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3, const double _u4, const double _u5);
// 	void ReadStartingPoints(const std::size_t i_sol, const double _u0, const double _u1, const double _u2, const double _u3);
// 	void ReadStartingPoints(const std::size_t i_sol, const double _u0, const double _u1);
// 	void ReadStartingPoints(const TArray2D u_in);
// 	void ReadGuess(const std::size_t i_sol, const double u_g);
// 	void UpdateStartingPoints(const TArray u);

// 	double GetExplicitPart(const std::size_t i_sol);
// 	double GetPhi(const std::size_t i_sol, const TArray u);

// 	void Step();

// 	void Step2(){ ; };
// 	void SetOrder(const std::size_t);
// 	std::size_t GetOrder();
// 	//size_t GetIterationsNumber(const size_t i_sol);
// };

#endif
