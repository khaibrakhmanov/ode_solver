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