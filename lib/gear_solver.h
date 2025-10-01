#pragma once
#include "ode_solver.h"

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