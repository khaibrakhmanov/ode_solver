#pragma once
#include "ctest.h"

class CTest1 : public CTest
{
public:
	CTest1();
	~CTest1();

protected:
	void SetEqs();
	void SetParams();
};