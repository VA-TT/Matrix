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
 o----\-o
 | \    |
 |   \  |
 o 0   \o 1  
///    ///
#endif

// Model Parameters
namespace modelParameters {
// Problem dimension: Considering the 2D implementation first
constexpr Index d{3};

// Unit vectors
Vector<double> i1{1.0, 0.0, 0.0};
Vector<double> i2{0.0, 1.0, 0.0};
Vector<double> i3{0.0, 0.0, 1.0};

// Geometry of the truss
constexpr Index nNodes{6}; // number of nodes
constexpr Index nBars{8};  // number of bars

Vector<Vector<double>> nodes{{0.0, 0.0, 0.0},  {10.0, 0.0, 0.0},
                             {0.0, 0.0, 10.0}, {10.0, 0.0, 10.0},
                             {0.0, 0.0, 20.0}, {10.0, 0.0, 20.0}};
// Bar connectivity: store node indices for each bar's origin and end
Vector<Index> barOrigin{0, 2, 1, 2, 2, 3, 3, 4};
Vector<Index> barEnd{2, 1, 3, 3, 4, 4, 5, 5};
Vector<Vector<double>> vectorBars(nBars), unitVectorBars(nBars);
Vector<double> lengthBars(nBars);
Vector<Index> nodeImposed{0, 1};
Vector<Index> nodeFree(nNodes - nodeImposed.size());

// Section dimension
double youngModulus{25e9}; // Module Young
double b{0.02}, h{0.05};
double A{b * h};
double alpha{youngModulus * A};

//  External loads at nodes
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
  assert(nodes.size() == nNodes && "numbers of nodes must be consistent!");
  assert(barOrigin.size() == nBars && barEnd.size() == nBars &&
         "numbers of bars must be consistent!");
  for (Index b{0}; b < nBars; ++b) {
    // compute bar vector from node coordinates
    vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
    lengthBars[b] = magnitude(vectorBars[b]);
    unitVectorBars[b] = vectorBars[b] / lengthBars[b];
  }
  std::cout << unitVectorBars;
  // Convert degrees to radians
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
  for (Index i = 0; i < forceExternal.size(); ++i) {
    f1[i] = forceExternal[i] * std::sin(thetaRadian[i]) * i1[0];
    f2[i] = -forceExternal[i] * std::cos(thetaRadian[i]) * i2[1];
    externalForce[i] = f1[i] + f2[i];
  }

  // Setting up
  int iteration{0};
  Vector<double> U(nNodes * d), UImposed(nNodes * d), UFree(nNodes * d);

  // Identify free nodes
  Vector<Index> nodeFree(nNodes - nodeImposed.size());
  Index nf{0};
  for (Index n{0}; n < nNodes; ++n) {
    bool isImposed = false;
    for (Index k{0}; k < nodeImposed.size(); ++k) {
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

  assert(nf == nodeFree.size() && "Mismatch in free node count");

  // Rigidity in small deformation configuration: N = k e (u2 - u1)
  Vector<double> k(nBars);
  for (Index b{0}; b < nBars; ++b) {
    k[b] = youngModulus * A / lengthBars[b];
  }

  std::cout << k;
  Vector<Matrix<double, d, d>> elementaryApplicationK(nBars);
  for (Index b{0}; b < nBars; ++b) {
    elementaryApplicationK[b] =
        k[b] * tensorProduct<d, d>(unitVectorBars[b], unitVectorBars[b]);
  }
  std::cout << elementaryApplicationK;

  // Matrix<double, 2, 2> connectivityMatrix{1.0, -1.0, -1.0, 1.0};
  // Vector<Matrix<double, d * 2, d * 2>> elementaryK(nBars);
  // for (Index b{0}; b < nBars; ++b) {
  //   elementaryK[b] =
  //       tensorProduct(connectivityMatrix, elementaryApplicationK[b]);
  // }
  // std::cout << elementaryK;

  // Connectivity Matrix
  Vector<Matrix<double, d, d * nNodes>> C(nNodes);

  const auto Id = Matrix<double, d, d>::identity();
  for (Index n{0}; n < nNodes; ++n) {
    Matrix<double, d, d * nNodes> Ci{}; // initialized to zero
    for (Index i = 0; i < d; ++i) {
      for (Index j = 0; j < d; ++j) {
        Ci(i, n * d + j) = Id(i, j);
      }
    }
    C[n] = Ci;
    // std::cout << C[n].getCols() << "  and  " << C[n].getRows() << std::endl;
  }

  Matrix<double, d * nNodes, d * nNodes> assemblyStiffnessK{};

  for (Index b{0}; b < nBars; ++b) {
    assemblyStiffnessK += (C[barEnd[b]] - C[barOrigin[b]]).transpose() *
                          elementaryApplicationK[b] *
                          (C[barEnd[b]] - C[barOrigin[b]]);
  }
  // std::cout << "Checkpoint: finished assembly loop" << std::endl;
  // // N += U; // This line causes a dimension mismatch error

  // return 0;
}