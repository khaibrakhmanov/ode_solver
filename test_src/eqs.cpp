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

double TExpODE::RHS(const double t, const double y) const
{
    return -GetParameter("alpha") * y;
}

double TExpODE::Solution(const double t) const
{
    return this->ic * std::exp(-GetParameter("alpha") * t);
}

// -------------------------------------------------------------

TLinODE::TLinODE(const double a, const double b, const double ic)
{
    this->ic = ic;
    this->SetParameter("a", a);
    this->SetParameter("b", b);
}

TLinODE::~TLinODE()
{
}

double TLinODE::RHS(const double t, const double y) const
{
    return GetParameter("a") -GetParameter("b") * y;
}

double TLinODE::Solution(const double t) const
{
    auto tau = 1 / GetParameter("b");

    return GetParameter("a") / GetParameter("b") * (1.0 - std::exp(-t / tau));
}

// -------------------------------------------------------------

TQuadODE::TQuadODE(const double a, const double b, const double ic)
{
    this->ic = ic;
    this->SetParameter("a", a);
    this->SetParameter("b", b);
}

TQuadODE::~TQuadODE()
{
}

double TQuadODE::RHS(const double t, const double y) const
{
    return GetParameter("a") - GetParameter("b") * std::pow(y, 2.0);
}

double TQuadODE::Solution(const double t) const
{
    auto tau = 0.5 / std::sqrt(GetParameter("a") * GetParameter("b"));

    return std::sqrt(GetParameter("a") / GetParameter("b")) * (1.0 - std::exp(-t / tau)) / (1.0 + std::exp(-t / tau));
}