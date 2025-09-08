#include "eqs.h"
#include <cmath>

TExpODE::TExpODE(const double alpha, const double ic)
{
    this->ic = ic;
    this->SetParameter("alpha", alpha);
}

TExpODE::~TExpODE()
{
}

double TExpODE::RHS(const double t, const double x) const
{
    return -GetParameter("alpha") * x;
}

double TExpODE::Solution(const double t, const double x) const
{
    return this->ic * std::exp(-GetParameter("alpha") * t);
}
