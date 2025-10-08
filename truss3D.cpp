#include "Matrix.h" //Approximative Comparsion
#include "Vector.h"
#include "dualDiffrentiation.h"
#include "vector"
#include <cassert> // for assert
#include <fstream> //working with files
#include <iomanip> //tab
#include <iostream>
#include <numbers>     // for std::numbers::pi
#include <stdexcept>   //throw exception
#include <type_traits> // precision

// Model
#if 0
 4      5
 o------o
 | \    |
 | 2 \  | 3
 o------o
 | \    |
 |   \  |
 o 0   \o 1  
///    ///
#endif

// Model Parameters
namespace modelParameters {
// Problem dimension: Considering the 2D implementation first
std::size_t d{2};

// Unit vectors
Vector<double> i1{1.0, 0.0};
Vector<double> i2{0.0, 1.0};

// Geometry of the truss
constexpr std::size_t nNodes{6}; // number of nodes
constexpr std::size_t nBars{8};  // number of bars

Vector<Vector<double>> N{{0.0, 0.0},   {10.0, 0.0}, {0.0, 10.0},
                         {10.0, 10.0}, {0, 20.0},   {10.0, 20.0}};
Vector<Vector<double>> barOrigin{N[0], N[2], N[1], N[2],
                                 N[2], N[3], N[3], N[4]};
Vector<Vector<double>> barEnd{N[2], N[1], N[3], N[3], N[4], N[4], N[5], N[5]};
Vector<Vector<double>> vectorBars(nBars), unitVectorBars(nBars);
Vector<double> lengthBars(nBars);
Vector<Index> nodeImposed{0, 1};
Vector<Index> nodeFree(nNodes - nodeImposed.size());

// Section dimension
double youngModulus{25e9}; // Module Young
double b{0.02}, h{0.05};
double A{b * h};
double alpha{youngModulus * A};

// External loads at nodes
Vector<Vector<double>> externalForce(nNodes);
Vector<double> forceExternal{0, 15, 0, 15, 15, 15}; // External Force
Vector<double> thetaDegree{90, 90, 90, 90,
                           90, 90}; // Inclined angle of the force with relative
                                    // to vertical converted to radian)

// Initialization of externalForce will be done in main()

} // namespace modelParameters

int main() {
  using namespace modelParameters;
  // Checking the input
  assert(N.size() == nNodes && "numbers of nodes must be consistent!");
  assert(barOrigin.size() == nBars && barEnd.size() == nBars &&
         "numbers of bars must be consistent!");
  for (Index b{0}; b < nBars; ++b) {
    vectorBars[b] = barEnd[b] - barOrigin[b];
    lengthBars[b] = magnitude(vectorBars[b]);
    unitVectorBars[b] = vectorBars[b] / lengthBars[b];
  }

  // Convert degrees → radians
  for (auto &angle : thetaDegree) {
    angle *= (std::numbers::pi / 180.0);
  }

  // Calculate external forces
  assert(forceExternal.size() == thetaDegree.size() &&
         thetaDegree.size() == nNodes &&
         "Numbers of nodal forces must respect number of nodes");
  Vector<double> thetaRadian{thetaDegree};
  Vector<double> f1(forceExternal.size()), f2(forceExternal.size()),
      externalForce(forceExternal.size());
  for (std::size_t i = 0; i < forceExternal.size(); ++i) {
    f1[i] = forceExternal[i] * std::sin(thetaRadian[i]) * i1[0];
    f2[i] = -forceExternal[i] * std::cos(thetaRadian[i]) * i2[1];
    externalForce[i] = f1[i] + f2[i];
  }

  // Setting up
  int iteration{0};
  Vector<double> U(nNodes * d), UImposed(nNodes * d), UFree(nNodes * d);

  // Xác định các nút tự do đúng cách
  Vector<Index> nodeFree(nNodes - nodeImposed.size());
  Index nf{0};
  for (Index n{0}; n < static_cast<Index>(nNodes); ++n) {
    bool isImposed = false;
    for (Index k{0}; k < static_cast<Index>(nodeImposed.size()); ++k) {
      if (n == nodeImposed[k]) {
        isImposed = true;
        break;
      }
    }
    if (!isImposed) {
      nodeFree[nf] = n;
      ++nf;
    }
  }

  assert(nf == static_cast<Index>(nodeFree.size()) &&
         "Mismatch in free node count");

    // N += U; // This line causes a dimension mismatch error

  std::cout << N;
  return 0;
}