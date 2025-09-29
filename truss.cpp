#include "Matrix.h" //Approximative Comparsion
#include "Vector.h"
#include "dualDiffrentiation.h"
#include "iostream"
#include "newtonRaphson.h"
#include <iomanip>     //tab
#include <stdexcept>   //throw exception
#include <type_traits> // precision
// Model
#if 0
|o----o-> F
     / 
    /
   /
|o/
#endif
// Model Parameters
namespace modelParameters {
double E{25e9}; // Module Youn
double l01{10}; // length
double l02{10}; // length
double F{15};   // Imposed Force
// Section dimension
double b1{0.2};
double h1{0.5};
double b2{0.2};
double h2{0.5};
double A1{b1 * h1};
double A2{b2 * h2};
double alpha1{E * A1};
double alpha1{E * A2};
// Unit vectors
Vector i1{1, 0};
Vector i2{0, 1};
} // namespace modelParameters

inline double constitutiveLaw(double alpha, double l, double l0) {
  return (alpha * (l - l0) / l0);
  return (alpha * (l * l - l0 * l0) / (2 * l0 * l0));
  return (alpha * std::log(l / l0));
}

int main() { return 0; }