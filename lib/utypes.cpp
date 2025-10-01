#include "utypes.h"

void AllocateVector1D(TArray *v, const size_t N)
{
	(*v).resize(N);
}

void AllocateVector2D(TArray2D *v, const size_t N1, const size_t N2) {
	(*v).resize(N1);
	for (size_t i = 0; i < N1; i++)
	{
		(*v)[i].resize(N2);
	};
};