#pragma once

#include <vector>
#include <map>
#include <string>

/// @brief class: Basic Ordinary differential equation
class TBasicODE
{
protected:
    /// @brief dictionary: ODE parameters (to be initialised)
    std::map <std::string, double> paramsDict;

public:
    TBasicODE();
    ~TBasicODE();

    /// @brief Set ODE parameter from the input
    /// @param name - parameter name, string
    /// @param val - parameter value, real
    void SetParameter(const std::string name, const double val);

    /// @brief Get ODE parameter value by name
    /// @param name - parameter name, string
    /// @return parameter value
    double GetParameter(const std::string name) const;
};

/// @brief class: single ODE
class TSingleODE : public TBasicODE
{
protected:

    /// @brief initial condition
    double ic;

public:
    TSingleODE(/* args */);
    ~TSingleODE();

    /// @brief Right-hand-side of the equation (value of the derivative)
    /// @param t - time (dependent variable)
    /// @param y - independent variable
    /// @return dy/dt
    virtual double RHS(const double t, const double y) const = 0;

    /// @brief Analytical solution of the ODE
    /// @param t - time (dependent variable)
    /// @return y(t)
    virtual double Solution(const double t) const = 0;

};

typedef double (*PFunction1)(double);
typedef double (*PFunction)(double, std::vector<double>);
//typedef TFunction* PFunction;
typedef std::vector<PFunction> TVectorFunction;
typedef std::vector<PFunction1> TVectorFunction1;

/// @brief class: system of ODEs
class TODEs : public TBasicODE
{
protected:

    std::vector<double> ic;
    TVectorFunction rhs;
    TVectorFunction1 solution;

public:
    TODEs(/* args */);
    ~TODEs();

    void AddRHS(const PFunction fun);
    void AddSolution(const PFunction1 sol);
    void AddIC(const double val);

    TVectorFunction* GetRHSPtr();
    TVectorFunction1* GetSolution();

};
