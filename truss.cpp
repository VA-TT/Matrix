#include "Matrix.h" //Approximative Comparsion
#include "Vector.h"
#include "dualDiffrentiation.h"
#include "vector"
#include <iomanip> //tab
#include <iostream>
#include <numbers>     // for std::numbers::pi
#include <stdexcept>   //throw exception
#include <type_traits> // precision

// Model
#if 0
|o-----o
      / |\
 ^i2 /  |  \> F
 |  /    theta
|o/---> i1
#endif

// Model Parameters
namespace modelParameters {
constexpr std::size_t n = 2; // number of equilibrium equations

// Unit vectors
Vector<double> i1{1.0, 0.0};
Vector<double> i2{0.0, 1.0};

// Size parameters
double E{25e9}; // Module Young
double a{10.0}; // Distance between node 0 and 1 = bar length

Vector<double> x0{(i1 * a + i2 * a)};
double l01{magnitude(x0)};          // length bar 1
double l02{magnitude(x0 - a * i2)}; // length bar 2
double force{15};                   // Imposed Force
double theta{20 * std::numbers::pi /
             180}; // Inclined angle of the force with relative to vertical

Vector<double> f1{force * std::sin(theta) * i1};  // horizontal component
Vector<double> f2{-force * std::cos(theta) * i2}; // vertical component
Vector<double> externalForce{f1 + f2};

// Section dimension
double b1{0.2};
double h1{0.5};
double b2{0.2};
double h2{0.5};
double A1{b1 * h1};
double A2{b2 * h2};
double alpha1{E * A1};
double alpha2{E * A2};

double epsilon{1e-6}; // tolerance
} // namespace modelParameters

// Original constitutive law
inline double constitutiveLaw(double alpha, double l, double l0) {
  return (alpha * (l - l0) / l0);
  // return (alpha * (l * l - l0 * l0) / (2 * l0 * l0));
  // return (alpha * std::log(l / l0));
}

// Template overload for Dual numbers
template <typename T>
inline auto constitutiveLaw(double alpha, T l, double l0) {
  return (alpha * (l - l0) / l0);
}

int main() {
  using namespace modelParameters;

  // Define func as a lambda that matches the expected signature for
  // automaticDiff
  auto func1 = [](auto x1) { return constitutiveLaw(alpha1, x1, l01); };
  auto func2 = [](auto x2) { return constitutiveLaw(alpha2, x2, l02); };

  // Setting up
  int iteration{0};
  int max_iteration{100};
  std::vector<int> iteration_array;
  std::vector<Vector<double>> deltaX_array;
  Vector<double> deltaX(2);
  Vector<double> F(2);

  std::cout << "=== Newton-Raphson Iteration ===" << std::endl;

  // Use Newton_Raphson method to find the approximation with tolerance
  while (iteration < max_iteration) {
    // Calculate updated informations
    Vector<double> x{x0 + deltaX};

    double l1{magnitude(x)};
    double l2{magnitude(x - a * i2)};
    Vector<double> e1{x / l1};
    Vector<double> e2{(x - a * i2) / l2};

    // Calculate F at iteration k
    F = -func1(l1) * e1 - func2(l2) * e2 + externalForce;
    std::cout << "Iteration " << iteration << ": |F| = " << magnitude(F)
              << std::endl;

    // Use F as the quality of approximation
    if (approximatelyEqualAbsRel(magnitude(F), epsilon))
      break;

    // Calculate NablaF at itaration k
    auto I = Matrix<double, 2, 2>::identity();
    Matrix<double, n, n> nablaF{};
    nablaF += -automaticDiff(func1, l1) * tensorProduct<n, n>(e1, e1);
    nablaF += -automaticDiff(func2, l2) * tensorProduct<n, n>(e2, e2);

    nablaF += -func1(l1) / l1 * (I - tensorProduct<n, n>(e1, e1));
    nablaF += -func2(l2) / l2 * (I - tensorProduct<n, n>(e2, e2));
    Vector<double> deltaX_increment = solveLinearSystem2(nablaF, -F);
    deltaX += deltaX_increment;

    // Saving the output
    iteration_array.push_back(iteration);
    deltaX_array.push_back(deltaX);

    iteration++;
  }

  if (iteration >= max_iteration) {
    std::cout << "Failed to converge after " << max_iteration << " iterations"
              << std::endl;
  }

  std::cout << "\nFinal displacement: " << deltaX << std::endl;
  std::cout << "Total iterations: " << iteration << std::endl;

  // Luu iteration va gia tri delta x
  return 0;
}