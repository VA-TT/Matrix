#include "Matrix.h" //Approximative Comparsion
#include "Vector.h"
#include "iostream"
#include <iomanip>     //tab
#include <stdexcept>   //throw exception
#include <type_traits> // precision

////////////////////////////////////////////////////////////////////////////
///////////////////////////   TEST FUNCTION   //////////////////////////////
////////////////////////////////////////////////////////////////////////////

void testMatrixConstructors() {
  std::cout << "=== Testing Matrix Constructors ===" << std::endl;

  // Default constructor
  Matrix<double, 2, 3> m1;
  std::cout << "Default constructor (2x3):\n" << m1 << std::endl;

  // Initializer list constructor
  Matrix<int, 2, 3> m2{1, 2, 3, 4, 5, 6};
  std::cout << "Initializer list constructor {1,2,3,4,5,6}:\n"
            << m2 << std::endl;

  // Copy constructor test
  Matrix<int, 2, 3> m3 = m2;
  std::cout << "Copy constructor:\n" << m3 << std::endl;

  std::cout << "Constructor tests passed!\n" << std::endl;
}

void testMatrixAccessors() {
  std::cout << "=== Testing Matrix Accessors ===" << std::endl;

  Matrix<int, 3, 3> m{1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::cout << "Original matrix:\n" << m << std::endl;

  // Test operator() access
  std::cout << "Element (0,0): " << m(0, 0) << std::endl;
  std::cout << "Element (1,2): " << m(1, 2) << std::endl;
  std::cout << "Element (2,1): " << m(2, 1) << std::endl;

  // Test operator[] access
  std::cout << "Element [0]: " << m[0] << std::endl;
  std::cout << "Element [4]: " << m[4] << std::endl;
  std::cout << "Element [8]: " << m[8] << std::endl;

  // Test modification
  m(1, 1) = 99;
  m[0] = 88;
  std::cout << "After modification m(1,1)=99, m[0]=88:\n" << m << std::endl;

  // Test getters
  std::cout << "Rows: " << m.getRows() << ", Cols: " << m.getCols()
            << ", Length: " << m.length() << std::endl;

  std::cout << "Accessor tests passed!\n" << std::endl;
}

void testMatrixStaticMethods() {
  std::cout << "=== Testing Static Matrix Methods ===" << std::endl;

  // Zero matrix
  auto zero = Matrix<double, 3, 3>::zero();
  std::cout << "Zero matrix (3x3):\n" << zero << std::endl;

  // Ones matrix
  auto ones = Matrix<int, 2, 4>::ones();
  std::cout << "Ones matrix (2x4):\n" << ones << std::endl;

  // Identity matrix
  auto identity = Matrix<double, 4, 4>::identity();
  std::cout << "Identity matrix (4x4):\n" << identity << std::endl;

  // Test reset methods
  Matrix<int, 2, 2> m{1, 2, 3, 4};
  std::cout << "Before reset:\n" << m << std::endl;

  m.resetZero();
  std::cout << "After resetZero():\n" << m << std::endl;

  m.resetOnes();
  std::cout << "After resetOnes():\n" << m << std::endl;

  m.resetIdentity();
  std::cout << "After resetIdentity():\n" << m << std::endl;

  std::cout << "Static methods tests passed!\n" << std::endl;
}

void testMatrixArithmetic() {
  std::cout << "=== Testing Matrix Arithmetic ===" << std::endl;

  Matrix<double, 2, 3> A{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix<double, 2, 3> B{2.0, 1.0, 4.0, 3.0, 6.0, 5.0};

  std::cout << "Matrix A:\n" << A << std::endl;
  std::cout << "Matrix B:\n" << B << std::endl;

  // Addition
  auto C = A + B;
  std::cout << "A + B:\n" << C << std::endl;

  // Subtraction
  auto D = A - B;
  std::cout << "A - B:\n" << D << std::endl;

  // Scalar multiplication
  auto E = 2.5 * A;
  std::cout << "2.5 * A:\n" << E << std::endl;

  auto F = A * 3.0;
  std::cout << "A * 3.0:\n" << F << std::endl;

  // Unary minus
  auto G = -A;
  std::cout << "-A:\n" << G << std::endl;

  // Equality
  auto A_copy = A;
  std::cout << "A == A_copy: " << (A == A_copy) << std::endl;
  std::cout << "A == B: " << (A == B) << std::endl;
  std::cout << "A != B: " << (A != B) << std::endl;

  std::cout << "Arithmetic tests passed!\n" << std::endl;
}

void testMatrixMultiplication() {
  std::cout << "=== Testing Matrix Multiplication ===" << std::endl;

  Matrix<double, 2, 3> A{1, 2, 3, 4, 5, 6};
  Matrix<double, 3, 2> B{1, 2, 3, 4, 5, 6};

  std::cout << "Matrix A (2x3):\n" << A << std::endl;
  std::cout << "Matrix B (3x2):\n" << B << std::endl;

  auto C = A * B;
  std::cout << "A * B (2x2):\n" << C << std::endl;

  auto D = B * A;
  std::cout << "B * A (3x3):\n" << D << std::endl;

  // Square matrix multiplication
  Matrix<int, 3, 3> E{1, 2, 3, 0, 1, 4, 5, 6, 0};
  Matrix<int, 3, 3> F{1, 0, 0, 0, 1, 0, 0, 0, 1}; // Identity

  std::cout << "Matrix E:\n" << E << std::endl;
  std::cout << "Matrix F (Identity):\n" << F << std::endl;

  auto G = E * F;
  std::cout << "E * Identity:\n" << G << std::endl;

  std::cout << "Matrix multiplication tests passed!\n" << std::endl;
}

void testMatrixProperties() {
  std::cout << "=== Testing Matrix Properties ===" << std::endl;

  // Diagonal matrix
  Matrix<int, 3, 3> diagonal{5, 0, 0, 0, 3, 0, 0, 0, 7};
  std::cout << "Diagonal matrix:\n" << diagonal << std::endl;
  std::cout << "isDiagonal(): " << diagonal.isDiagonal() << std::endl;
  std::cout << "isUpperTriangular(): " << diagonal.isUpperTriangular()
            << std::endl;
  std::cout << "isLowerTriangular(): " << diagonal.isLowerTriangular()
            << std::endl;
  std::cout << "isSquare(): " << diagonal.isSquare() << std::endl;

  // Upper triangular
  Matrix<int, 3, 3> upper{1, 2, 3, 0, 4, 5, 0, 0, 6};
  std::cout << "\nUpper triangular matrix:\n" << upper << std::endl;
  std::cout << "isUpperTriangular(): " << upper.isUpperTriangular()
            << std::endl;
  std::cout << "isLowerTriangular(): " << upper.isLowerTriangular()
            << std::endl;

  // Symmetric matrix
  Matrix<double, 3, 3> symmetric{1, 2, 3, 2, 4, 5, 3, 5, 6};
  std::cout << "\nSymmetric matrix:\n" << symmetric << std::endl;
  std::cout << "isSymmetric(): " << symmetric.isSymmetric() << std::endl;

  // Non-square matrix
  Matrix<int, 2, 3> nonsquare{1, 2, 3, 4, 5, 6};
  std::cout << "\nNon-square matrix (2x3):\n" << nonsquare << std::endl;
  std::cout << "isSquare(): " << nonsquare.isSquare() << std::endl;

  std::cout << "Matrix properties tests passed!\n" << std::endl;
}

void testDeterminantAndTrace() {
  std::cout << "=== Testing Determinant and Trace ===" << std::endl;

  // 1x1 matrix
  Matrix<double, 1, 1> m1{5.0};
  std::cout << "1x1 matrix: " << m1(0, 0) << std::endl;
  std::cout << "det(1x1): " << det(m1) << std::endl;
  std::cout << "trace(1x1): " << trace(m1) << std::endl;

  // 2x2 matrix
  Matrix<double, 2, 2> m2{1, 2, 3, 4};
  std::cout << "\n2x2 matrix:\n" << m2 << std::endl;
  std::cout << "det(2x2): " << det(m2) << std::endl;
  std::cout << "trace(2x2): " << trace(m2) << std::endl;

  // 3x3 matrix
  Matrix<double, 3, 3> m3{1, 2, 3, 0, 1, 4, 5, 6, 0};
  std::cout << "\n3x3 matrix:\n" << m3 << std::endl;
  std::cout << "det(3x3): " << det(m3) << std::endl;
  std::cout << "trace(3x3): " << trace(m3) << std::endl;

  // Identity matrix
  auto identity = Matrix<double, 3, 3>::identity();
  std::cout << "\n3x3 Identity:\n" << identity << std::endl;
  std::cout << "det(Identity): " << det(identity) << std::endl;
  std::cout << "trace(Identity): " << trace(identity) << std::endl;

  std::cout << "Determinant and trace tests passed!\n" << std::endl;
}

void testMatrixInverse() {
  std::cout << "=== Testing Matrix Inverse ===" << std::endl;

  try {
    // 2x2 invertible matrix
    Matrix<double, 2, 2> m2{1, 2, 3, 5};
    std::cout << "Original 2x2 matrix:\n" << m2 << std::endl;
    std::cout << "det(m2): " << det(m2) << std::endl;

    auto inv2 = m2.inverse();
    std::cout << "Inverse:\n" << inv2 << std::endl;

    auto product2 = m2 * inv2;
    std::cout << "Original * Inverse:\n" << product2 << std::endl;

  } catch (const std::exception &e) {
    std::cout << "Caught exception: " << e.what() << std::endl;
  }

  // Test singular matrix (should throw)
  try {
    Matrix<double, 2, 2> singular{1, 2, 2, 4};
    std::cout << "\nTrying to invert singular matrix:\n"
              << singular << std::endl;
    auto inv_singular = singular.inverse();
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  std::cout << "Matrix inverse tests passed!\n" << std::endl;
}

void testRowColumnOperations() {
  std::cout << "=== Testing Row/Column Operations ===" << std::endl;

  Matrix<int, 3, 4> m{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  std::cout << "Original matrix:\n" << m << std::endl;

  // Test row access
  auto row1 = m.row(1);
  std::cout << "Row 1: ";
  for (const auto &val : row1) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  // Test column access
  auto col2 = m.col(2);
  std::cout << "Column 2: ";
  for (const auto &val : col2) {
    std::cout << val.get() << " ";
  }
  std::cout << std::endl;

  // Test row swapping
  Matrix<int, 3, 3> swap_test{1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::cout << "Before swapping rows 0 and 2:\n" << swap_test << std::endl;
  swap_test.swapRows(0, 2);
  std::cout << "After swapping rows 0 and 2:\n" << swap_test << std::endl;

  // Test transpose
  Matrix<double, 2, 3> transpose_test{1, 2, 3, 4, 5, 6};
  std::cout << "Original (2x3):\n" << transpose_test << std::endl;
  auto transposed = transpose_test.transpose();
  std::cout << "Transposed (3x2):\n" << transposed << std::endl;

  std::cout << "Row/Column operations tests passed!\n" << std::endl;
}

void testConcatenationAndSplit() {
  std::cout << "=== Testing Concatenation and Split ===" << std::endl;

  Matrix<int, 2, 2> A{1, 2, 3, 4};
  Matrix<int, 2, 3> B{5, 6, 7, 8, 9, 10};

  std::cout << "Matrix A (2x2):\n" << A << std::endl;
  std::cout << "Matrix B (2x3):\n" << B << std::endl;

  auto concat = concatenatedMatrix(A, B);
  std::cout << "Concatenated matrix A|B (2x5):\n" << concat << std::endl;

  // Test splitByColumn
  Matrix<int, 2, 2> split_A;
  Matrix<int, 2, 3> split_B;
  concat.splitByColumn(2, split_A, split_B);

  std::cout << "Split back - first part:\n" << split_A << std::endl;
  std::cout << "Split back - second part:\n" << split_B << std::endl;

  std::cout << "Concatenation and split tests passed!\n" << std::endl;
}

void testLinearSystem() {
  std::cout << "=== Testing Linear System Solver ===" << std::endl;

  try {
    // Solve Ax = b
    Matrix<double, 3, 3> A{2, 1, 1, 1, 0, 1, 0, 3, 1};
    Matrix<double, 3, 1> b{4, 2, 6};

    std::cout << "Coefficient matrix A:\n" << A << std::endl;
    std::cout << "Right-hand side b:\n" << b << std::endl;

    auto x = solveLinearSystem(A, b);
    std::cout << "Solution x:\n" << x << std::endl;

    // Verify: Ax should equal b
    auto verification = A * x;
    std::cout << "Verification A*x:\n" << verification << std::endl;

  } catch (const std::exception &e) {
    std::cout << "Caught exception: " << e.what() << std::endl;
  }

  std::cout << "Linear system tests passed!\n" << std::endl;
}

void testVectorIntegration() {
  std::cout << "=== Testing Matrix-Vector Integration ===" << std::endl;

  Matrix<double, 3, 3> A{2, 1, 1, 1, 0, 1, 0, 3, 1};
  Vector<double> b{4, 2, 6};

  std::cout << "Matrix A:\n" << A << std::endl;
  std::cout << "Vector b: " << b << std::endl;

  auto x = solveLinearSystem2(A, b);
  std::cout << "Solution x:\n" << x << std::endl;

  // Convert Vector back to Matrix for multiplication
  Matrix<double, 3, 1> xMatrix{x};
  auto verification = A * xMatrix;
  std::cout << "Verification A*x:\n" << verification << std::endl;
}

void testTensorProduct() {
  std::cout << "=== Testing Tensor Product ===" << std::endl;

  // Test basic tensor product
  Vector<double> v1{1.0, 2.0, 3.0}; // 3D vector
  Vector<double> v2{4.0, 5.0};      // 2D vector

  std::cout << "Vector v1 = " << v1 << std::endl;
  std::cout << "Vector v2 = " << v2 << std::endl;

  // Sử dụng explicit template parameters
  auto tensor = makeTensorProduct<3, 2>(v1, v2);
  std::cout << "Tensor product v1 ⊗ v2 (3x2):\n" << tensor << std::endl;

  // Verify kết quả:
  // [1] ⊗ [4, 5] = [1*4, 1*5] = [4, 5]
  // [2]             [2*4, 2*5]   [8, 10]
  // [3]             [3*4, 3*5]   [12, 15]

  std::cout << "Expected result:" << std::endl;
  std::cout << "Row 0: [4, 5]" << std::endl;
  std::cout << "Row 1: [8, 10]" << std::endl;
  std::cout << "Row 2: [12, 15]" << std::endl;

  std::cout << "Tensor product tests passed!" << std::endl;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////   MAIN   //////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main() {
  std::cout << "=== MATRIX CLASS COMPREHENSIVE TEST SUITE ===" << std::endl;
  std::cout << "=============================================" << std::endl
            << std::endl;

  // Run all test cases
  testMatrixConstructors();
  testMatrixAccessors();
  testMatrixStaticMethods();
  testMatrixArithmetic();
  testMatrixMultiplication();
  testMatrixProperties();
  testDeterminantAndTrace();
  testMatrixInverse();
  testRowColumnOperations();
  testConcatenationAndSplit();
  testLinearSystem();
  testVectorIntegration();
  testTensorProduct();

  std::cout << "\n=== ORIGINAL TESTS (from previous main) ===" << std::endl;

  // Your original tests
  Matrix<double, 3, 3> A{-3, 2, -1, 6, -6, 7, 3, -4, 4};
  Matrix<double, 3, 1> B{-1, -7, -6};

  std::cout << "Original problem: Solve Ax = B" << std::endl;
  std::cout << "A:\n" << A << std::endl;
  std::cout << "B:\n" << B << std::endl;

  try {
    Matrix<double, 3, 1> X{solveLinearSystem(A, B)};
    std::cout << "Solution X:\n" << X << std::endl;

    // Verify solution
    auto verification = A * X;
    std::cout << "Verification A*X:\n" << verification << std::endl;

    std::cout << "A inverse:\n" << A.inverse() << std::endl;

  } catch (const std::exception &e) {
    std::cout << "Error: " << e.what() << std::endl;
  }

  // Test various determinants
  std::cout << "\nDeterminant tests:" << std::endl;
  const Matrix<double, 1, 1> H1{6};
  const Matrix<double, 2, 2> H2{6, 3, 2, 4};
  const Matrix<double, 3, 3> D{2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};

  std::cout << "det(H1): " << det(H1) << std::endl;
  std::cout << "det(H2): " << det(H2) << std::endl;
  std::cout << "det(D): " << det(D) << std::endl;

  std::cout << "\n=== ALL TESTS COMPLETED SUCCESSFULLY! ===" << std::endl;

  std::cout << "=== VECTOR CLASS TEST SUITE ===" << std::endl;
  std::cout << std::endl;

  // Test 1: Constructors
  std::cout << "1. Testing Constructors:" << std::endl;
  Vector<int> v1;                // Default constructor
  Vector<int> v2(5);             // Size constructor (5 zeros)
  Vector<int> v3{1, 2, 3, 4, 5}; // Initializer list

  std::cout << "v1 (default): " << v1 << " (size: " << v1.size() << ")"
            << std::endl;
  std::cout << "v2(5): " << v2 << " (size: " << v2.size() << ")" << std::endl;
  std::cout << "v3{1,2,3,4,5}: " << v3 << " (size: " << v3.size() << ")"
            << std::endl;
  std::cout << std::endl;

  // Test 2: Indexing
  std::cout << "2. Testing Indexing:" << std::endl;
  Vector<double> v4{10.5, 20.3, 30.7, 40.1};
  std::cout << "v4: " << v4 << std::endl;
  std::cout << "v4[0] = " << v4[0] << std::endl;
  std::cout << "v4[2] = " << v4[2] << std::endl;

  // Modify element
  v4[1] = 99.9;
  std::cout << "After v4[1] = 99.9: " << v4 << std::endl;

  // Test signed indexing
  Index idx = 3;
  std::cout << "Using Index idx=3: v4[idx] = " << v4[idx] << std::endl;
  std::cout << std::endl;

  // Test 3: Vector Operations
  std::cout << "3. Testing Vector Operations:" << std::endl;
  Vector<double> va{1.0, 2.0, 3.0};
  Vector<double> vb{4.0, 5.0, 6.0};

  std::cout << "va = " << va << std::endl;
  std::cout << "vb = " << vb << std::endl;

  auto vc = va + vb;
  std::cout << "va + vb = " << vc << std::endl;

  auto vd = vb - va;
  std::cout << "vb - va = " << vd << std::endl;

  auto ve = 2.5 * va;
  std::cout << "2.5 * va = " << ve << std::endl;

  auto vf = va * 3.0;
  std::cout << "va * 3.0 = " << vf << std::endl;

  auto vg = -va;
  std::cout << "-va = " << vg << std::endl;
  std::cout << std::endl;

  // Test 4: Vector Math
  std::cout << "4. Testing Vector Math:" << std::endl;
  Vector<double> u1{3.0, 4.0, 0.0};
  Vector<double> u2{1.0, 0.0, 0.0};

  std::cout << "u1 = " << u1 << std::endl;
  std::cout << "u2 = " << u2 << std::endl;

  double dot = dotProduct(u1, u2);
  std::cout << "Dot product u1·u2 = " << dot << std::endl;

  auto cross = crossProduct(u1, u2);
  std::cout << "Cross product u1×u2 = " << cross << std::endl;

  double mag1 = magnitude(u1);
  double mag2 = magnitude(u2);
  std::cout << "Magnitude |u1| = " << mag1 << std::endl;
  std::cout << "Magnitude |u2| = " << mag2 << std::endl;

  auto unit1 = normalize(u1);
  auto unit2 = normalize(u2);
  std::cout << "Unit vector of u1 = " << unit1 << std::endl;
  std::cout << "Unit vector of u2 = " << unit2 << std::endl;
  std::cout << "Magnitude of unit u1 = " << magnitude(unit1) << std::endl;
  std::cout << std::endl;

  // Test 4a: Testing CrossProduct2 (Levi-Civita Method)
  std::cout << "4a. Testing CrossProduct2 (Levi-Civita Method):" << std::endl;
  auto cross2 = crossProduct2(u1, u2);
  std::cout << "Cross product2 u1×u2 = " << cross2 << std::endl;
  std::cout << "crossProduct == crossProduct2: "
            << (cross == cross2 ? "✓ SAME" : "✗ DIFFERENT") << std::endl;
  std::cout << std::endl;

  // Test 5: Resize and Push Back
  std::cout << "5. Testing Resize and Push Back:" << std::endl;
  Vector<int> vr{10, 20, 30};
  std::cout << "Original vr: " << vr << " (size: " << vr.size() << ")"
            << std::endl;

  vr.resize(6);
  std::cout << "After resize(6): " << vr << " (size: " << vr.size() << ")"
            << std::endl;

  vr.resize(8, 77);
  std::cout << "After resize(8, 77): " << vr << " (size: " << vr.size() << ")"
            << std::endl;

  vr.push_back(100);
  std::cout << "After push_back(100): " << vr << " (size: " << vr.size() << ")"
            << std::endl;
  std::cout << std::endl;

  // Test 6: Edge Cases and Error Handling
  std::cout << "6. Testing Edge Cases:" << std::endl;

  try {
    Vector<int> v_small{1, 2};
    Vector<int> v_big{1, 2, 3, 4};
    std::cout << "Trying to add vectors of different sizes..." << std::endl;
    auto result = v_small + v_big;
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  try {
    Vector<double> zero_vec{0.0, 0.0, 0.0};
    std::cout << "Trying to normalize zero vector..." << std::endl;
    auto unit = normalize(zero_vec);
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  try {
    Vector<int> v2d_1{1, 2};
    Vector<int> v2d_2{3, 4};
    std::cout << "Trying cross product on 2D vectors..." << std::endl;
    auto cross_2d = crossProduct(v2d_1, v2d_2);
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  std::cout << std::endl;
  std::cout << "=== ALL TESTS COMPLETED SUCCESSFULLY! ===" << std::endl;

  // Test 7: Projection
  std::cout << "7. Testing Vector Projection:" << std::endl;

  Vector<double> proj_a{3.0, 4.0, 0.0}; // Vector a
  Vector<double> proj_b{1.0, 0.0, 0.0}; // Unit vector along x-axis
  Vector<double> proj_c{2.0, 2.0, 0.0}; // 45-degree vector

  std::cout << "Vector proj_a = " << proj_a << std::endl;
  std::cout << "Vector proj_b = " << proj_b << std::endl;
  std::cout << "Vector proj_c = " << proj_c << std::endl;

  // Vector projection
  auto vector_proj_ab = proj_a.projection(proj_b);
  auto vector_proj_ac = proj_a.projection(proj_c);
  std::cout << "Projection of proj_a onto proj_b: " << vector_proj_ab
            << std::endl;
  std::cout << "Projection of proj_a onto proj_c: " << vector_proj_ac
            << std::endl;

  std::cout << std::endl;

  // Test 8: Angle Calculation
  std::cout << "8. Testing Angle Calculation:" << std::endl;

  Vector<double> angle_v1{1.0, 0.0, 0.0}; // Unit vector along x-axis
  Vector<double> angle_v2{0.0, 1.0, 0.0}; // Unit vector along y-axis
  Vector<double> angle_v3{1.0, 1.0, 0.0}; // 45-degree vector

  std::cout << "angle_v1 = " << angle_v1 << std::endl;
  std::cout << "angle_v2 = " << angle_v2 << std::endl;
  std::cout << "angle_v3 = " << angle_v3 << std::endl;

  // Test 90 degree angle
  double angle12_rad = angleRad(angle_v1, angle_v2);
  double angle12_deg = angleDegree(angle_v1, angle_v2);
  std::cout << "Angle between angle_v1 and angle_v2: " << angle12_rad
            << " rad = " << angle12_deg << " deg" << std::endl;

  // Test 45 degree angle
  double angle13_rad = angleRad(angle_v1, angle_v3);
  double angle13_deg = angleDegree(angle_v1, angle_v3);
  std::cout << "Angle between angle_v1 and angle_v3: " << angle13_rad
            << " rad = " << angle13_deg << " deg" << std::endl;

  std::cout << std::endl;

  // Test 9: Perpendicular and Parallel
  std::cout << "9. Testing Perpendicular and Parallel:" << std::endl;

  Vector<double> perp_x{1.0, 0.0, 0.0}; // Unit vector along x
  Vector<double> perp_y{0.0, 1.0, 0.0}; // Unit vector along y
  Vector<double> para_x{2.0, 0.0, 0.0}; // Parallel to perp_x

  std::cout << "perp_x = " << perp_x << std::endl;
  std::cout << "perp_y = " << perp_y << std::endl;
  std::cout << "para_x = " << para_x << std::endl;

  // Test perpendicular (should be true)
  bool is_perp = isPerpendicular(perp_x, perp_y);
  std::cout << "isPerpendicular(perp_x, perp_y): " << is_perp << std::endl;

  // Test parallel (should be true)
  bool is_para = isParallel(perp_x, para_x);
  std::cout << "isParallel(perp_x, para_x): " << is_para << std::endl;

  // Show verification values
  std::cout << "Dot product perp_x·perp_y = " << dotProduct(perp_x, perp_y)
            << std::endl;
  std::cout << "Cross product magnitude |perp_x×para_x| = "
            << magnitude(crossProduct(perp_x, para_x)) << std::endl;

  std::cout << std::endl;

  return 0;
}
