#pragma once
#include "ctest.h"

class CTest2 : public CTest
{
public:
	CTest2();
	~CTest2();

protected:
	void SetEqs();
	void SetParams();
};