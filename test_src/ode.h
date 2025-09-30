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

/// @brief Type: pointer to a function of one real argument
typedef double (*PFunction)(double);
/// @brief vector of pointers to functions each depending on one real argument
typedef std::vector<PFunction> TVectorFunction;

/// @brief Type: pointer to a function of one real and one vector arguments
typedef double (*PFunctionVectArg)(double, std::vector<double>);
/// @brief Type: vector of pointers to functions each depending on real and one vector arguments
typedef std::vector<PFunctionVectArg> TVectorFunctionVectArg;

/// @brief class: system of ODEs
/// The system has a following form:
///   dy/dt = rhs,
/// where 
///   y   = {y0, y1, y2, ...} - vector of unknown functions (solution of the system)
///   rhs = {f0, f1, f2, ...} - vector-function of the derivatives in the system (a.k.a. Right-Hand-Side)
/// The system is closed by a vecot of initial conditions:
///   ic  = {ic0, ic1. ...}
class TODEs : public TBasicODE
{
protected:
    /// @brief vector of inital conditions
    std::vector<double> ic;
    /// @brief vector of "right-hand-side" functions (derivatives in the system)
    TVectorFunctionVectArg rhs;
    /// @brief vector-function of the exact solution of the system
    TVectorFunction solution;

public:
    TODEs(/* args */);
    ~TODEs();

    /// @brief Add RHS (derivative) to a vector of system equations
    /// @param fun function to be added
    void AddRHS(const PFunctionVectArg fun);
    /// @brief Add exact solution to a vector of solutions
    /// @param sol function to be added
    void AddSolution(const PFunction sol);
    /// @brief Add a value of the initial condition to a vector of initial conditions
    /// @param val initial condition to be added
    void AddIC(const double val);

    /// @brief Get pointer to a vector-function of the system equations
    /// @return *rhs
    TVectorFunctionVectArg* GetRHSPtr();
    /// @brief Get pointer to a vector of functions of the exact solutions
    /// @return *solution
    TVectorFunction* GetSolution();

};
