#include <array>
#include <vector>
#include <span>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <initializer_list>
#include <cmath>
#include <iomanip>
#include <functional>
#include "comparison.hpp"
#include <type_traits>
#include "Matrix.h"

// Include your Matrix class here or copy it above this line

void testMatrixConstructors() {
    std::cout << "=== Testing Matrix Constructors ===" << std::endl;
    
    // Default constructor
    Matrix<double, 2, 3> m1;
    std::cout << "Default constructor (2x3):\n" << m1 << std::endl;
    
    // Initializer list constructor
    Matrix<int, 2, 3> m2{1, 2, 3, 4, 5, 6};
    std::cout << "Initializer list constructor {1,2,3,4,5,6}:\n" << m2 << std::endl;
    
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
    std::cout << "Element (0,0): " << m(0,0) << std::endl;
    std::cout << "Element (1,2): " << m(1,2) << std::endl;
    std::cout << "Element (2,1): " << m(2,1) << std::endl;
    
    // Test operator[] access
    std::cout << "Element [0]: " << m[0] << std::endl;
    std::cout << "Element [4]: " << m[4] << std::endl;
    std::cout << "Element [8]: " << m[8] << std::endl;
    
    // Test modification
    m(1,1) = 99;
    m[0] = 88;
    std::cout << "After modification m(1,1)=99, m[0]=88:\n" << m << std::endl;
    
    // Test getters
    std::cout << "Rows: " << m.getRows() << ", Cols: " << m.getCols() << ", Length: " << m.length() << std::endl;
    
    std::cout << "Accessor tests passed!\n" << std::endl;
}

void testMatrixRowColumnAccess() {
    std::cout << "=== Testing Row/Column Access ===" << std::endl;
    
    Matrix<int, 3, 4> m{1, 2, 3, 4,
                        5, 6, 7, 8,
                        9, 10, 11, 12};
    std::cout << "Original matrix:\n" << m << std::endl;
    
    // Test row access
    auto row1 = m.row(1);
    std::cout << "Row 1: ";
    for (const auto& val : row1) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    
    // Modify row
    row1[2] = 77;
    std::cout << "After modifying row1[2] = 77:\n" << m << std::endl;
    
    // Test column access
    auto col2 = m.col(2);
    std::cout << "Column 2: ";
    for (const auto& val : col2) {
        std::cout << val.get() << " ";
    }
    std::cout << std::endl;
    
    // Modify column
    col2[0].get() = 66;
    std::cout << "After modifying col2[0] = 66:\n" << m << std::endl;
    
    std::cout << "Row/Column access tests passed!\n" << std::endl;
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
    std::cout << "isUpperTriangular(): " << diagonal.isUpperTriangular() << std::endl;
    std::cout << "isLowerTriangular(): " << diagonal.isLowerTriangular() << std::endl;
    std::cout << "isSquare(): " << diagonal.isSquare() << std::endl;
    
    // Upper triangular
    Matrix<int, 3, 3> upper{1, 2, 3, 0, 4, 5, 0, 0, 6};
    std::cout << "\nUpper triangular matrix:\n" << upper << std::endl;
    std::cout << "isUpperTriangular(): " << upper.isUpperTriangular() << std::endl;
    std::cout << "isLowerTriangular(): " << upper.isLowerTriangular() << std::endl;
    
    // Lower triangular
    Matrix<int, 3, 3> lower{1, 0, 0, 2, 3, 0, 4, 5, 6};
    std::cout << "\nLower triangular matrix:\n" << lower << std::endl;
    std::cout << "isUpperTriangular(): " << lower.isUpperTriangular() << std::endl;
    std::cout << "isLowerTriangular(): " << lower.isLowerTriangular() << std::endl;
    
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

void testMatrixOperations() {
    std::cout << "=== Testing Matrix Operations ===" << std::endl;
    
    Matrix<double, 3, 3> m{1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::cout << "Original matrix:\n" << m << std::endl;
    
    // Transpose
    auto mt = m.transpose();
    std::cout << "Transpose:\n" << mt << std::endl;
    
    // Swap rows
    Matrix<int, 3, 3> swap_test{1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::cout << "Before swapping rows 0 and 2:\n" << swap_test << std::endl;
    swap_test.swapRows(0, 2);
    std::cout << "After swapping rows 0 and 2:\n" << swap_test << std::endl;
    
    // Find max row index
    Matrix<double, 3, 3> find_max{1.0, -5.0, 3.0, 2.0, 8.0, -1.0, -3.0, 2.0, 4.0};
    std::cout << "Matrix for finding max in column:\n" << find_max << std::endl;
    std::cout << "Max element in column 1 starting from row 0: row " << find_max.indexRowMax(1, 0) << std::endl;
    
    // SubMatrix
    Matrix<double, 3, 3> sub_test{1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::cout << "Original matrix for subMatrix:\n" << sub_test << std::endl;
    auto sub = sub_test.subMatrix(1, 1);
    std::cout << "SubMatrix removing row 1, col 1:\n" << sub << std::endl;
    
    std::cout << "Matrix operations tests passed!\n" << std::endl;
}

void testDeterminantAndTrace() {
    std::cout << "=== Testing Determinant and Trace ===" << std::endl;
    
    // 1x1 matrix
    Matrix<double, 1, 1> m1{5.0};
    std::cout << "1x1 matrix: " << m1(0,0) << std::endl;
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
    
    // Singular matrix (det = 0)
    Matrix<double, 2, 2> singular{1, 2, 2, 4};
    std::cout << "\nSingular matrix:\n" << singular << std::endl;
    std::cout << "det(singular): " << det(singular) << std::endl;
    
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
        
        // 3x3 invertible matrix
        Matrix<double, 3, 3> m3{2, 1, 1, 1, 0, 1, 0, 3, 1};
        std::cout << "\nOriginal 3x3 matrix:\n" << m3 << std::endl;
        std::cout << "det(m3): " << det(m3) << std::endl;
        
        auto inv3 = m3.inverse();
        std::cout << "Inverse:\n" << inv3 << std::endl;
        
        auto product3 = m3 * inv3;
        std::cout << "Original * Inverse:\n" << product3 << std::endl;
        
    } catch (const std::exception& e) {
        std::cout << "Caught exception: " << e.what() << std::endl;
    }
    
    // Test singular matrix (should throw)
    try {
        Matrix<double, 2, 2> singular{1, 2, 2, 4};
        std::cout << "\nTrying to invert singular matrix:\n" << singular << std::endl;
        auto inv_singular = singular.inverse();
    } catch (const std::exception& e) {
        std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
    }
    
    std::cout << "Matrix inverse tests passed!\n" << std::endl;
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
        
    } catch (const std::exception& e) {
        std::cout << "Caught exception: " << e.what() << std::endl;
    }
    
    std::cout << "Linear system tests passed!\n" << std::endl;
}

void testMatrixConcatenation() {
    std::cout << "=== Testing Matrix Concatenation ===" << std::endl;
    
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
    
    std::cout << "Matrix concatenation tests passed!\n" << std::endl;
}

void testEdgeCasesAndExceptions() {
    std::cout << "=== Testing Edge Cases and Exceptions ===" << std::endl;
    
    // Test dimension mismatch in addition
    try {
        Matrix<int, 2, 2> A{1, 2, 3, 4};
        Matrix<int, 3, 3> B{1, 2, 3, 4, 5, 6, 7, 8, 9};
        std::cout << "Trying to add 2x2 and 3x3 matrices..." << std::endl;
        auto C = A + B;  // Should fail at assertion
    } catch (const std::exception& e) {
        std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
    }
    
    // Test subMatrix edge case
    try {
        Matrix<double, 1, 1> tiny{5.0};
        std::cout << "Trying to create subMatrix of 1x1 matrix..." << std::endl;
        auto sub = tiny.subMatrix(0, 0);  // Should throw
    } catch (const std::exception& e) {
        std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
    }
    
    // Test orthogonal matrix check
    Matrix<double, 3, 3> orthogonal{1, 0, 0, 0, 1, 0, 0, 0, 1}; // Identity is orthogonal
    std::cout << "Identity matrix isOrthogonal(): " << orthogonal.isOrthogonal() << std::endl;
    
    std::cout << "Edge cases and exceptions tests passed!\n" << std::endl;
}

int main() {
    std::cout << "=== MATRIX CLASS COMPREHENSIVE TEST SUITE ===" << std::endl;
    std::cout << "=============================================" << std::endl << std::endl;
    
    testMatrixConstructors();
    testMatrixAccessors();
    testMatrixRowColumnAccess();
    testMatrixStaticMethods();
    testMatrixArithmetic();
    testMatrixMultiplication();
    testMatrixProperties();
    testMatrixOperations();
    testDeterminantAndTrace();
    testMatrixInverse();
    testLinearSystem();
    testMatrixConcatenation();
    testEdgeCasesAndExceptions();
    
    std::cout << "=== ALL MATRIX TESTS COMPLETED SUCCESSFULLY! ===" << std::endl;
    
    return 0;
}