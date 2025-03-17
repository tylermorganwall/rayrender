#ifndef ASSERTH
#define ASSERTH

#ifndef NDEBUG
#define ASSERT(condition) (assertCondition((condition), #condition))
#else
#define ASSERT(condition) ((void)0)
#endif

// Function to throw error if condition is false
void assertCondition(bool condition, const char* conditionStr);

#endif