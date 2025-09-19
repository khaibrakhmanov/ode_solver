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

/// @brief Pointer to a function of one real argument
typedef double (*PFunction)(double);
/// @brief vector of pointers to functions each depending on one real argument
typedef std::vector<PFunction> TVectorFunction;

/// @brief Pointer to a function of one real and one vector arguments
typedef double (*PFunctionVectArg)(double, std::vector<double>);
/// @brief vector of pointers to functions each depending on real and one vector arguments
typedef std::vector<PFunctionVectArg> TVectorFunctionVectArg;

/// @brief class: system of ODEs
class TODEs : public TBasicODE
{
protected:

    std::vector<double> ic;
    TVectorFunctionVectArg rhs;
    TVectorFunction solution;

public:
    TODEs(/* args */);
    ~TODEs();

    void AddRHS(const PFunctionVectArg fun);
    void AddSolution(const PFunction sol);
    void AddIC(const double val);

    TVectorFunctionVectArg* GetRHSPtr();
    TVectorFunction* GetSolution();

};
