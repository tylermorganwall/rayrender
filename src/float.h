#ifndef FLOATH
#define FLOATH

// #define RAY_FLOAT_AS_DOUBLE

#ifdef RAY_FLOAT_AS_DOUBLE
typedef double Float;
#else
typedef float Float;
#endif 

#endif