#pragma once
#include "ode_solver.h"

class TBasicRungeKuttaSolver : public TBaseODESolver
{
protected:
	/// @brief Matrix of the RK coefficients
	TArray2D coeff;

	/// intermediate solution: y1 = y0 + 0.5 * C0
	TArray y1;
	/// intermediate solution: y2 = y1 + 0.5 * C1
	TArray y2;
	/// intermediate solution: y3 = y2 + 0.5 * C2
	TArray y3;
	/// @brief Calculate intermediate vectors at current time step
	/// @param y0 vector of initial values
	void CalcIntemediateVariables(const TArray y0);

public:
	TBasicRungeKuttaSolver();

	/// @brief Calculate the coefficients of the scheme
	/// @param y0 vector of the initial values
	/// @param _dt current time step
	/// @param _time surrent time
	virtual void CalcCoefficients(const TArray y0, const double _dt, const double _time) = 0;

	//virtual void Step() = 0; //with h
	//virtual void Step2() = 0;// with h/2 + h/2

	/// @brief Returns the RK coefficient for given indeces
	/// @param i index of the coefficient
	/// @param j index of the equation
	/// @return coeff[i, j]
	double GetC(const std::size_t i, const std::size_t j);

	/// @brief Set the RK coefficient for given indeces
	/// @param i index of the coefficient
	/// @param j index of the equation
	/// @param arg value of the coefficient
	void SetC(const std::size_t i, const std::size_t j, const double arg);

	/// <summary>
	/// Returns pointer to the array of the scheme coefficients
	/// </summary>
	/// <returns>pointer to coeff</returns>
	TArray2D* GetCoefficientsPtr();

	/// @brief Allocate memory for the scheme coefficients
	void AllocateCoefficients();

	/// @brief Get the order of the scheme
	/// @return order
	std::size_t GetOrder();

	/// @brief Get the order of the scheme
	/// @param o order
	void SetOrder(const std::size_t o);
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
