#include "ode.h"

TBasicODE::TBasicODE()
{

}

TBasicODE::~TBasicODE()
{

}

void TBasicODE::SetParameter(const std::string name, const double val)
{
    paramsDict.insert({name, val});
}

double TBasicODE::GetParameter(const std::string name) const
{
    return paramsDict.at(name);
}

TSingleODE::TSingleODE()
{
    ic = 0.0;
}

TSingleODE::~TSingleODE()
{
    
}

TODEs::TODEs()
{
}

TODEs::~TODEs()
{
}

void TODEs::AddRHS(const PFunction fun)
{
    rhs.push_back(fun);
    //rhs.at(i) = fun;
}

void TODEs::AddSolution(const PFunction1 sol)
{
    solution.push_back(sol);
    //solution.at(i) = fun;
}

void TODEs::AddIC(const double val)
{
    ic.push_back(val);
}

TVectorFunction* TODEs::GetRHSPtr()
{
    return &rhs;
}

TVectorFunction* TODEs::GetSolution()
{
    return &solution;
}
