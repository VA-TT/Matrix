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
#include <functional> // std::reference_wrapper
#include "comparison.hpp"

template <typename T, std::size_t nRows, std::size_t nCols>
class Matrix
{
private:
  std::array<T, nRows * nCols> m_elements{};

public:
  // Constructors, Destructors
  Matrix(std::initializer_list<T> list)
  {
    assert(list.size() == nRows * nCols);
    std::copy(list.begin(), list.end(), m_elements.begin());
  }
  Matrix() = default;
  Matrix(const Matrix &) = default;
  Matrix(Matrix &&) = default;
  Matrix &operator=(const Matrix &) = default;
  Matrix &operator=(Matrix &&) = default;
  ~Matrix() = default;

  // 0 Matrix, 1 Matrix and Identity Matrix and functions to reset, static to make them create one copy only
  static Matrix zero()
  {
    return Matrix{};
  }
  void resetZero()
  {
    (*this) = Matrix::zero();
  }

  static Matrix ones()
  {
    Matrix result{};
    std::fill(result.m_elements.begin(), result.m_elements.end(), 1);
    return result;
  }
  void resetOnes()
  {
    (*this) = Matrix::ones();
  }

  static Matrix identity()
  {
    Matrix result{};
    for (std::size_t i = 0; i < std::min(nRows, nCols); ++i)
      result(i, i) = 1.0;
    return result;
  }
  void resetIdentity()
  {
    (*this) = Matrix::identity();
  }

  // Accessing the elements in the array with one parameter (i)
  T &operator[](std::size_t i)
  {
    assert(i < this->length());
    return m_elements[i];
  }
  const T &operator[](std::size_t i) const
  {
    assert(i < this->length());
    return m_elements[i];
  }

  // Accessing the elements in the matrices with two parameters (i,j)
  T &operator()(std::size_t i, std::size_t j)
  {
    assert(i < nRows && j < nCols);
    return m_elements[i * nCols + j];
    // Consider return m_elements[(i - 1) * nCols + (j - 1)]; in order to accessing the matrices with index starting from 1 in mathematic
  }

  const T &operator()(std::size_t i, std::size_t j) const
  {
    assert(i < nRows && j < nCols);
    return m_elements[i * nCols + j];
  }

  // Getters to get number of rows & columns + total elements numbers
  constexpr std::size_t getCols() const { return nCols; }
  constexpr std::size_t getRows() const { return nRows; }
  constexpr int length() const { return nRows * nCols; }

  // Const Reference to a whole row(i) or column(j)
  // std::array<T, nCols> row(std::size_t i) const
  // {
  //   assert(i < nRows);
  //   std::array<T, nCols> row_i{};
  //   for (std::size_t j{0}; j < nCols; ++j)
  //     row_i[j] = (*this)(i, j);
  //   return row_i;
  // }

  // std::array<T, nRows> col(std::size_t j) const
  // {
  //   assert(j < nCols);
  //   std::array<T, nRows> col_j{};
  //   for (std::size_t i{0}; i < nRows; ++i)
  //     col_j[i] = (*this)(i, j);
  //   return col_j;
  // }

  // Reference to row(i)
  std::span<T> row(std::size_t i)
  {
    assert(i < nRows);
    return std::span<T>(&m_elements[i * nCols], nCols);
  }

  // Reference to col(j): couldn't use std::span (data is not continuous in m_element) and std::array here (as with array reference_wrapper default constructor will fail)
  std::vector<std::reference_wrapper<T>> col(std::size_t j)
  {
    assert(j < nCols);
    std::vector<std::reference_wrapper<T>> col_refs;
    col_refs.reserve(nRows);
    for (std::size_t i = 0; i < nRows; ++i)
    {
      col_refs.push_back(std::ref((*this)(i, j)));
    }
    return col_refs;
  }

  // outputting matrix
  friend std::ostream &operator<<(std::ostream &out, const Matrix &matrix)
  {
    constexpr int tab = 6;
    for (std::size_t i = 0; i < nRows; ++i)
    {
      out << "|";
      for (std::size_t j = 0; j < nCols; ++j)
      {
        out << std::setw(tab) << matrix(i, j);
      }
      out << " |" << '\n';
    }
    return out;
  }

  // Matrix algebra operators
  friend Matrix operator*(T k, const Matrix &m)
  {
    Matrix<T, nRows, nCols> result{m};
    for (auto &e : result.m_elements)
    {
      e *= k;
    }
    return result;
  }
  friend Matrix operator*(const Matrix &m, T k)
  {
    return k * m;
  }

  friend Matrix operator+(const Matrix &m1, const Matrix &m2)
  {
    assert(m1.getCols() == m2.getCols() && m1.getRows() == m2.getRows() && "Unable to perform matrix addition/substraction.");
    Matrix<T, nRows, nCols> result{m1};
    for (std::size_t i{0}; i < result.length(); ++i)
      result[i] += m2[i];
    return result;
  }

  friend Matrix operator-(const Matrix &m1, const Matrix &m2)
  {
    return m1 + (-1) * m2;
  }

  Matrix operator-() const
  {
    return (-1) * (*this);
  }

  friend bool operator==(const Matrix &m1, const Matrix &m2) // co the vut ra ngoai duoc
  {
    double tolerance = 5e-4; // Should make it relatively
    if (m1.getCols() != m2.getCols() || m1.getRows() != m2.getRows())
    {
      return false;
    }
    // for (std::size_t i{0}; i < m1.length(); i++)
    // {
    //   if (std::abs(m1[i] - m2[i]) > tolerance)
    //   {
    //     return false;
    //   }
    // }
    // return true;
    for (std::size_t i{0}; i < m1.length(); i++)
    {
      if (!approximatelyEqualAbsRel(m1[i], m2[i]))
      {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const Matrix &m1, const Matrix &m2)
  {
    return !(m1 == m2);
  }

  // Separate a matrix to two matrices by a column (position at the beginning of the second matrix)
  template <std::size_t leftCols, std::size_t rightCols>
  void splitByColumn(std::size_t colPos,
                     Matrix<T, nRows, leftCols> &A,
                     Matrix<T, nRows, rightCols> &B) const
  {
    assert(leftCols == colPos);
    assert(rightCols == nCols - colPos);

    for (std::size_t i = 0; i < nRows; ++i)
      for (std::size_t j = 0; j < nCols; ++j)
      {
        if (j < colPos)
          A(i, j) = (*this)(i, j);
        else
          B(i, j - colPos) = (*this)(i, j);
      }
  }

  // // Pointer to a whole row(i) or column(j)
  // std::array<T *, nCols> rowPointer(std::size_t i)
  // {
  //   assert(i < nRows);
  //   std::array<T *, nCols> row_pointer{};
  //   for (std::size_t j = 0; j < nCols; ++j)
  //   {
  //     row_pointer[j] = &(*this)(i, j);
  //   }
  //   return row_pointer;
  // }

  // std::array<T *, nRows> colPointer(std::size_t j)
  // {
  //   assert(j < nCols);
  //   std::array<T *, nRows> col_pointer{};
  //   for (std::size_t i = 0; i < nRows; ++i)
  //   {
  //     col_pointer[i] = &(*this)(i, j);
  //   }
  //   return col_pointer;
  // }

  // Swap 2 rows
  void swapRows(std::size_t i1, std::size_t i2)
  {
    for (std::size_t j{0}; j < nCols; ++j)
    {
      std::swap((*this)(i1, j), (*this)(i2, j));
    }
  }

  // Find the row with the highest number of elements that aren't 0
  int MostZeroRow(double tolerance = 5e-4)
  {
    std::array<int, nRows> exceed0{};
    std::size_t index{0};
    for (std::size_t i{0}; i < nRows; i++)
    {
      for (const auto &element : (*this).row(i))
      {
        if (std::abs(element) < tolerance)
          ++exceed0[i];
      }
    }
    for (std::size_t i{1}; i < nRows; i++)
    {
      if (exceed0[i] > exceed0[index])
      {
        index = i;
      }
    }
    // Tim index chua gia tri max trong array exceed0
    return static_cast<int>(index);
  }

  // Find the row with the max element at a given column (Noted below pivot)
  std::size_t indexRowMax(std::size_t col, std::size_t startRow = 0) const
  {
    assert(col < nCols && startRow < nRows);
    std::size_t maxRow = startRow;
    T maxValue = std::abs((*this)(startRow, col));
    for (std::size_t i = startRow + 1; i < nRows; ++i)
    {
      T value = std::abs((*this)(i, col));
      if (value > maxValue)
      {
        maxValue = value;
        maxRow = i;
      }
    }
    return maxRow;
  }

  // Transpose matrix
  Matrix<T, nCols, nRows> transpose() const
  {
    Matrix<T, nCols, nRows> result{};
    for (std::size_t j = 0; j < nCols; ++j)
    {
      for (std::size_t i = 0; i < nRows; ++i)
      {
        result(j, i) = (*this)(i, j);
      }
    }
    return result;
  }

  // Inverse Matrix
  Matrix<T, nRows, nCols> inverse() const
  {
    if (!isSquare())
      throw std::invalid_argument("Inverse only defined for square matrices.");
    // static_assert(nRows == nCols, "Inverse only defined for square matrices.");
    // int maxCount{100};
    // int count{0};
    // bool computeFlag = false;
    Matrix<T, nRows, nCols> inverseMatrix{};
    Matrix<T, nRows, nCols> testIdentity{};
    Matrix<T, nRows, (nCols + nCols)> augmentedMatrix{concatenatedMatrix(*this, Matrix<T, nRows, nCols>::identity())};
    // while (count < maxCount && !computeFlag)
    for (std::size_t pivotIndex{0}; pivotIndex < nRows; ++pivotIndex)
    {
      std::cout << "Pivot index " << pivotIndex << std::endl;
      std::size_t maxIndex{augmentedMatrix.indexRowMax(pivotIndex, pivotIndex)};
      if (maxIndex != pivotIndex)
      {

        augmentedMatrix.swapRows(pivotIndex, maxIndex);
        std::cout << "Swap rows " << pivotIndex << " and " << maxIndex << std::endl;
        std::cout << augmentedMatrix << std::endl;
      }
      T pivot{augmentedMatrix(pivotIndex, pivotIndex)};
      if (approximatelyEqualAbsRel(pivot, 0.0)) // approximate to 0
      {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
      }
      if (!approximatelyEqualAbsRel(pivot, 1.0))
      {
        augmentedMatrix.row(pivotIndex) = augmentedMatrix.row(pivotIndex) / pivot; // normalize pivot
      }
      for (std::size_t i{0}; i < nRows; ++i)
      {
        if (i == pivotIndex)
          continue;
        // If the element is already 0.0, move on
        if (!approximatelyEqualAbsRel(augmentedMatrix(i, pivotIndex), 0.0))
        {
          T factor = -augmentedMatrix(i, pivotIndex); // /pivot = /1.0
          augmentedMatrix.row(i) = augmentedMatrix.row(i) + factor * augmentedMatrix.row(pivotIndex);
          std::cout << "Multiply row " << pivotIndex << " by " << factor << " and add it to row " << i << std::endl;
          std::cout << augmentedMatrix << std::endl;
        }
      }
    }
    augmentedMatrix.splitByColumn(nCols, testIdentity, inverseMatrix);
    if (testIdentity == Matrix::identity())
      return inverseMatrix;
    else
      throw std::runtime_error("Matrix inversion failed: result is not an identity matrix.");
  }

  // Bool function to get the characteristics of matrix
  bool isDiagonal() const
  {
    for (std::size_t i = 0; i < nRows; ++i)
    {
      for (std::size_t j = 0; j < nCols; ++j)
      {
        if (i != j)
        {
          if ((*this)(i, j) != 0)
          {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool isUpperTriangular() const
  {
    for (std::size_t i = 0; i < nRows; ++i)
    {
      for (std::size_t j = 0; j < nCols; ++j)
      {
        if (i > j)
        {
          if ((*this)(i, j) != 0)
          {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool isLowerTriangular() const
  {
    for (std::size_t i = 0; i < nRows; ++i)
    {
      for (std::size_t j = 0; j < nCols; ++j)
      {
        if (i < j)
        {
          if ((*this)(i, j) != 0)
          {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool isSquare() const
  {
    return nRows == nCols;
  }

  bool isSymmetric() const
  {
    return ((*this) == this->transpose());
  }

  bool isSkewSymmetric() const
  {
    return (-(*this) == (this->transpose()));
  }

  bool isOrthogonal() const
  {
    if (this->isSquare())
    {
      return (this->inverse() == this->transpose());
    }
    return false;
  }

  // Resize function?
};

//////////////////////////////////   END OF MATRIX CLASS   //////////////////////////////////

// Multiply a scalar to a row
// template <typename T, std::size_t nCols>
// std::array<T, nCols> operator*(const T &k, const std::array<T, nCols> &row)
// {
//   std::array<T, nCols> result{};
//   for (std::size_t i = 0; i < nCols; ++i)
//     result[i] = row[i] * k;
//   return result;
// }

// template <typename T, std::size_t nCols>
// std::array<T, nCols> operator*(const std::array<T, nCols> &row, const T &k)
// {
//   return k * row;
// }

// template <typename T, std::size_t nCols>
// std::array<T, nCols> operator/(const std::array<T, nCols> &row, const T &k)
// {
//   std::array<T, nCols> result{};
//   for (std::size_t i = 0; i < nCols; ++i)
//     result[i] = row[i] / k;
//   return result;
// }

// template <typename T, std::size_t nCols>
// std::array<T, nCols> operator/(const T &k, const std::array<T, nCols> &row)
// {
//   std::array<T, nCols> result{};
//   for (std::size_t i = 0; i < nCols; ++i)
//     result[i] = k / row[i];
//   return result;
// }
// Change the value of row directly
template <typename T>
std::span<T> operator*(std::span<T> row, const T &k)
{
  for (auto &element : row)
    element *= k;
  return row;
}

template <typename T>
std::span<T> operator*(const T &k, std::span<T> row)
{

  return row * k;
}

template <typename T>
std::span<T> operator/(std::span<T> row, const T &k)
{
  for (auto &element : row)
    element /= k;
  return row;
}

template <typename T>
std::span<T> operator/(const T &k, std::span<T> row)
{
  for (auto &element : row)
    element = k / element;
  return row;
}

template <typename T, std::size_t nCols>
std::span<T> operator+(std::span<T> row1, std::span<T> row2)
{
  for (std::size_t j{0}; j < nCols; j++)
  {
    row1[j] += row2[j];
  }
  return row1;
}

template <typename T, std::size_t nCols>
std::span<T> operator-(std::span<T> row1, std::span<T> row2)
{
  for (std::size_t j{0}; j < nCols; j++)
  {
    row1[j] -= row2[j];
  }
  return row1;
}

// Adding a row multiplied with a const to another row
// template <typename T, std::size_t nCols>
// const std::array<T, nCols> operator+(const std::array<T, nCols> &row1, const std::array<T, nCols> &row2)
// {
//   std::array<T, nCols> result{};
//   for (std::size_t j{0}; j < nCols; j++)
//   {
//     result[j] = row1[j] + row2[j];
//   }
//   return result;
// }

template <typename T, std::size_t R1, std::size_t C1, std::size_t R2, std::size_t C2>
Matrix<T, R1, (C1 + C2)> concatenatedMatrix(const Matrix<T, R1, C1> &A, const Matrix<T, R2, C2> &B)
{
  assert(R1 == R2 && "The number of rows of both matrices must match to create an concatenated matrix.\n");
  Matrix<T, R1, (C1 + C2)> concatenatedMatrix{};
  for (std::size_t i = 0; i < R1; ++i)
    for (std::size_t j = 0; j < (C1 + C2); ++j)
    {
      if (j < C1)
      {
        assert(j < A.getCols() && "Error index exceeding matrix's size.\n");
        concatenatedMatrix(i, j) = A(i, j);
      }
      else
      {
        assert((j - C1) < B.getCols() && "Error index exceeding matrix's size.\n");
        concatenatedMatrix(i, j) = B(i, j - C1);
      }
    }
  return concatenatedMatrix;
}

// Gaussian elimination to solve Ax = B
template <typename T, std::size_t R1, std::size_t C1>
Matrix<T, R1, 1> gaussianElimination(const Matrix<T, R1, C1> &A, const Matrix<T, R1, 1> &B)
{
  // assert(R1 == R2 && "Not suitable for solving Ax = B.");
  constexpr std::size_t augCols{C1 + 1};
  Matrix<T, R1, (augCols)> augmentedMatrix{concatenatedMatrix(A, B)};
  T coefficient{};

  for (std::size_t i = 0; i < R1 - 1; ++i)
  {
    T pivot{augmentedMatrix(i, i)};
    if (pivot == 0)
    {
      std::cout << "Row interchange must first be perfomed.\n";
    }
    for (std::size_t k{i + 1}; k < R1; ++k)
    {
      coefficient = augmentedMatrix(k, i) / pivot;
      for (std::size_t j = 0; j < augCols; ++j)
      {
        augmentedMatrix(k, j) -= (coefficient * augmentedMatrix(i, j));
      }
    }
  }
  // Gaussian elimination to solve for x in Ax = B
  Matrix<T, R1, 1> result{};
  std::cout << augmentedMatrix << '\n';
  return result;
}
// Reduced Row echilon form to calculate matrix inverse (to be implemented)

// LU Decomposition (to be implemented)

// Matrix multiplication
template <typename T, std::size_t R1, std::size_t C1, std::size_t R2, std::size_t C2>
Matrix<T, R1, C2> operator*(const Matrix<T, R1, C1> &m1, const Matrix<T, R2, C2> &m2)
{
  assert(C1 == R2 && "Not suitable for matrix multiplication.");
  Matrix<T, R1, C2> result{};
  for (std::size_t i = 0; i < R1; ++i)
    for (std::size_t j = 0; j < C2; ++j)
      for (std::size_t k = 0; k < C1; ++k)
        result(i, j) += m1(i, k) * m2(k, j);
  return result;
}

// Orthogonal bool function
template <typename T, std::size_t R1, std::size_t C1, std::size_t R2, std::size_t C2>
bool arePairOrthogonal(const Matrix<T, R1, C1> &m1, const Matrix<T, R2, C2> &m2)
{
  if ((m1 != Matrix<T, R1, C1>::zero()) && (m2 != Matrix<T, R2, C2>::zero()))
  {
    return ((m1.transpose() * m2) == Matrix<T, R1, C2>::identity());
  }
  return false; // Both matrices must not be 0
}

// Calculate trace(A)
template <typename T, std::size_t nRows, std::size_t nCols>
T trace(const Matrix<T, nRows, nCols> &m)
{
  static_assert(nRows == nCols, "Trace applies only to squared matrix");
  T result{};
  for (std::size_t i = 0; i < nRows; ++i)
    result += m(i, i);
  return result;
}

//////////////////////////////////   MAIN   //////////////////////////////////
int main()
{

  // Test: Trace cho ma trận vuông
  Matrix<int, 3, 3> C{1, 2, 3,
                      4, 5, 6,
                      7, 8, 9};
  std::cout << "Trace(C): " << trace(C) << std::endl;

  Matrix<double, 3, 3> D{2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};
  std::cout << D.inverse() << '\n';
  return 0;
}