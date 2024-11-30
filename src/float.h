#ifndef FLOATH
#define FLOATH

// #define RAY_FLOAT_AS_DOUBLE

#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 


#include <cmath>

template<class T>
inline T ffmin(T a, T b) { return(a < b ? a : b);}

template<class T>
inline T ffmax(T a, T b) { return(a > b ? a : b);}


template<>
inline float ffmin(float a, float b) { return(fmin(a,b));}

template<class T>
inline double ffmax(double a, double b) { return(fminf(a,b));}


template<class T>
inline T ffabs(T a) { return(a > 0 ? a : -a);}

template<>
inline float ffabs(float a) { return(fabsf(a));}

template<>
inline double ffabs(double a) { return(fabs(a));}


#endif