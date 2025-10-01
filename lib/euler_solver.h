#pragma once
#include "ode_solver.h"

/// @brief Class: basic Euler solver
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
