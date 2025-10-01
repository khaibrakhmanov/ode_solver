#pragma once
#include <vector>

using namespace std;

/// @brief user type: one-dimensional vector of real values
typedef std::vector<double> TArray;

/// @brief @brief user type: two-dimensional vector of real values
typedef std::vector<TArray> TArray2D;

/// @brief Allocate memory for a 1-dimensional array
/// @param v array to allocate memory for
/// @param N number of elements
void AllocateVector1D(TArray *v, const size_t N);

/// @brief Allocate memory for a 2-dimensional array
/// @param v array to allocate memory for
/// @param N1 number of elements in the 1-st dimension
/// @param N2 number of elements in the 2-nd dimension
void AllocateVector2D(TArray2D *v, const size_t N1, const size_t N2);
