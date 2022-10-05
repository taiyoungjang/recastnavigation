/**
@defgroup detour Detour

Members in this module are wrappers around the standard math library
*/

#ifndef DETOURMATH_H
#define DETOURMATH_H

#include <math.h>
// This include is required because libstdc++ has problems with isfinite
// if cmath is included before math.h.
#include <cmath>

inline double dtMathFabsf(double x) { return fabs(x); }
inline double dtMathSqrtf(double x) { return sqrt(x); }
inline double dtMathFloorf(double x) { return floor(x); }
inline double dtMathCeilf(double x) { return ceil(x); }
inline double dtMathCosf(double x) { return cos(x); }
inline double dtMathSinf(double x) { return sin(x); }
inline double dtMathAtan2f(double y, double x) { return atan2(y, x); }
inline bool dtMathIsfinite(double x) { return std::isfinite(x); }

#endif
