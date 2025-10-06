#include "Matrix.h"
#include "Vector.h"
#include "dualDiffrentiation.h"
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace modelParameters {
constexpr std::size_t n = 2;

Vector<double> i1{1.0, 0.0};
Vector<double> i2{0.0, 1.0};

double E{25e9};
double a{10.0};

Vector<double> x0{i1 * a + i2 * a};
double l01{magnitude(x0)};
double l02{magnitude(x0 - a * i2)};
double force{15};
double theta{20 * std::numbers::pi / 180.0};

Vector<double> f1{force * std::sin(theta) * i1};
Vector<double> f2{-force * std::cos(theta) * i2};
Vector<double> externalForce{f1 + f2};

double b1{0.2}, h1{0.5};
double b2{0.2}, h2{0.5};
double A1{b1 * h1}, A2{b2 * h2};
double alpha1{E * A1}, alpha2{E * A2};

double epsilon{1e-6};
int max_iteration{100};
} // namespace modelParameters

inline double constitutiveLaw(double alpha, double l, double l0) {
  return alpha * (l - l0) / l0;
}
template <typename T>
inline auto constitutiveLaw(double alpha, T l, double l0) {
  return alpha * (l - l0) / l0;
}

// Nếu chưa có (giữ nguyên nếu đã tồn tại ở nơi khác)
inline Vector<double> solveLinearSystem2(const Matrix<double, 2, 2> &A,
                                         const Vector<double> &b) {
  double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
  if (std::abs(det) < 1e-14)
    throw std::runtime_error("Singular matrix");
  Vector<double> x(2);
  x[0] = (b[0] * A(1, 1) - b[1] * A(0, 1)) / det;
  x[1] = (A(0, 0) * b[1] - A(1, 0) * b[0]) / det;
  return x;
}

int main() {
  using namespace modelParameters;

  auto func1 = [](auto l) { return constitutiveLaw(alpha1, l, l01); };
  auto func2 = [](auto l) { return constitutiveLaw(alpha2, l, l02); };

  Vector<double> deltaX(2);
  Vector<double> F(2);

  std::cout << "=== Newton-Raphson Iteration ===\n";

  int iteration = 0;
  while (iteration < max_iteration) {
    Vector<double> x = x0 + deltaX;
    double l1 = magnitude(x);
    double l2 = magnitude(x - a * i2);

    if (l1 < 1e-14 || l2 < 1e-14)
      throw std::runtime_error(
          "Bar length too small (possible singular config)");

    Vector<double> e1 = x / l1;
    Vector<double> e2 = (x - a * i2) / l2;

    F = -func1(l1) * e1 - func2(l2) * e2 + externalForce;
    double Fn = magnitude(F);
    std::cout << "Iter " << iteration << " |F| = " << Fn << "\n";

    if (Fn < epsilon) {
      std::cout << "Converged.\n";
      break;
    }

    auto I = Matrix<double, 2, 2>::identity();
    Matrix<double, n, n> nablaF{};

    // Debug derivatives
    auto dT1_auto = automaticDiff(func1, l1);
    auto dT2_auto = automaticDiff(func2, l2);

    // Fallback analytic derivatives (for linear constitutive law)
    double dT1 = std::abs(dT1_auto) < 1e-30 ? (alpha1 / l01) : dT1_auto;
    double dT2 = std::abs(dT2_auto) < 1e-30 ? (alpha2 / l02) : dT2_auto;

    // Assemble (using fallback derivatives)
    nablaF += -dT1 * tensorProduct<n, n>(e1, e1);
    nablaF += -dT2 * tensorProduct<n, n>(e2, e2);
    double T1 = func1(l1);
    double T2 = func2(l2);
    nablaF += -(T1 / l1) * (I - tensorProduct<n, n>(e1, e1));
    nablaF += -(T2 / l2) * (I - tensorProduct<n, n>(e2, e2));

    // Compute determinant for debug
    double detJ = nablaF(0, 0) * nablaF(1, 1) - nablaF(0, 1) * nablaF(1, 0);

    std::cout << "  l1=" << l1 << " l2=" << l2 << "\n";
    std::cout << "  e1=(" << e1[0] << "," << e1[1] << ") "
              << " e2=(" << e2[0] << "," << e2[1] << ")\n";
    std::cout << "  T1=" << T1 << " T2=" << T2 << " dT1=" << dT1_auto
              << " dT2=" << dT2_auto << " [used dT1=" << dT1 << " dT2=" << dT2
              << "]\n";
    std::cout << "  nablaF=\n" << nablaF << "  det(nablaF)=" << detJ << "\n";

    if (std::abs(detJ) < 1e-20) {
      throw std::runtime_error("Singular Jacobian (det≈0) after assembly");
    }

    Vector<double> incr = solveLinearSystem2(nablaF, -F);
    deltaX = deltaX + incr;

    iteration++;
  }

  if (iteration >= max_iteration)
    std::cout << "Failed to converge in " << max_iteration << " iterations\n";

  std::cout << "Final displacement: " << deltaX << "\n";
  std::cout << "Final position: " << x0 + deltaX << "\n";
  return 0;
}