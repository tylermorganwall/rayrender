#include "assert.h"
#include <string>
#include <stdexcept>

// Function to throw error if condition is false
void assertCondition(bool condition, const char* conditionStr) {
  if (!condition) {
    std::string message = "Assertion failed: ";
    message += conditionStr;
    throw std::runtime_error(message);
  }
}