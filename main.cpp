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
#include <type_traits> // thêm dòng này

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
    std::ios oldState(nullptr);
    oldState.copyfmt(out); // save stream state

    if constexpr (std::is_floating_point_v<T>)
    {
      out << std::fixed << std::setprecision(3);
    }

    constexpr int tab = 10;
    for (std::size_t i = 0; i < nRows; ++i)
    {
      out << "|";
      for (std::size_t j = 0; j < nCols; ++j)
      {
        out << std::setw(tab) << matrix(i, j);
      }
      out << " |" << '\n';
    }

    out.copyfmt(oldState);
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
    Matrix<T, nRows, nCols> inverseMatrix{};
    Matrix<T, nRows, nCols> testIdentity{};
    Matrix<T, nRows, (nCols + nCols)> augmentedMatrix{concatenatedMatrix(*this, Matrix<T, nRows, nCols>::identity())};

    for (std::size_t pivotIndex{0}; pivotIndex < nRows; ++pivotIndex)
    {
      // std::cout << "Pivot index " << pivotIndex << std::endl;
      std::size_t maxIndex{augmentedMatrix.indexRowMax(pivotIndex, pivotIndex)};
      if (maxIndex != pivotIndex)
      {
        augmentedMatrix.swapRows(pivotIndex, maxIndex);
        // std::cout << "Swap rows " << pivotIndex << " and " << maxIndex << std::endl;
        // std::cout << augmentedMatrix << std::endl;
      }
      T pivot{augmentedMatrix(pivotIndex, pivotIndex)};
      if (approximatelyEqualAbsRel(pivot, 0.0))
      {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
      }
      auto pivotSpan = augmentedMatrix.row(pivotIndex);
      // Normalize pivot in place
      if (!approximatelyEqualAbsRel(pivot, 1.0))
      {
        for (auto &val : pivotSpan)
          val /= pivot;
        // std::cout << "Normalize row " << pivotIndex << " by pivot = " << pivot << "\n"
        //           << augmentedMatrix << '\n';
      }

      // Eliminate other rows
      for (std::size_t i{0}; i < nRows; ++i)
      {
        // Skip the pivot row
        if (i == pivotIndex)
          continue;

        T factor = augmentedMatrix(i, pivotIndex); // no need to /pivot as pivot = 1.0 already
        if (approximatelyEqualAbsRel(factor, 0.0)) // Skip value under or above the pilot that are 0.0
          continue;

        // std::cout << "Eliminate row " << i << " using pivot row " << pivotIndex
        //           << " (factor = " << factor << ")\n";
        auto targetSpan = augmentedMatrix.row(i);
        for (std::size_t j = 0; j < pivotSpan.size(); ++j)
        {
          targetSpan[j] -= factor * pivotSpan[j];
        }
        // std::cout << augmentedMatrix << '\n';
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

// Solve Ax = B by taking x = A.inverse()*B
template <typename T, std::size_t R1, std::size_t C1>
Matrix<T, R1, 1> solveLinearSystem(const Matrix<T, R1, C1> &A, const Matrix<T, R1, 1> &B)
{
  try
  {
    Matrix<T, R1, C1> inverseA{A.inverse()};
    Matrix<T, R1, 1> result{inverseA * B};
    return result;
  }
  catch (const std::exception &)
  {
    throw std::invalid_argument("A is not invertible.");
  }
}

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
  Matrix<double, 3, 3> A{-3, 2, -1, 6, -6, 7, 3, -4, 4};
  Matrix<double, 3, 1> B{-1, -7, -6};
  Matrix<double, 3, 1> X{solveLinearSystem(A, B)};
  std::cout << "X = " << '\n'
            << X << '\n';
  Matrix<double, 3, 3> C{-3, 2, -1,
                         6, -6, 7,
                         3, -4, 4};
  std::cout << C.inverse() << '\n';

  const Matrix<double, 3, 3> D{2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};
  std::cout << D.inverse() << '\n';
  const Matrix<double, 5, 5> E{2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 9.0, 5.0, 3.0, 2.0, 6.0, 2.0, 4.0, 6.0, 5.0, 1.0, 1.0, 7.0, 5.0, 2.0, 3.0};
  std::cout << E.inverse() << '\n';
  Matrix<double, 5, 5> F{-7.0, -3.0, -7.0, -6.0, -1.0, 3.0, 0.0, 2.0, -4.0, -3.0, -9.0, -10.0, -2.0, 9.0, -3.0, -2.0, 7.0, -5.0, -10.0, -10.0, -1.0, 1.0, 0.0, 8.0, -5.0};
  std::cout << F.inverse() << '\n';
  Matrix<double, 5, 5> G{-7.06711357635508, 2.86204331759976, 8.72551512419021, -3.49596266237113, 1.09637006931382,
                         2.29334230427082, -8.86162756725844, -5.53715539554703, -7.73341321854147, 8.0623899548081,
                         -7.06854274385278, -8.29455270479789, 8.44531043955609, 4.84975158392336, -8.9356105292099,
                         4.54902139094058, 6.89377456364728, 3.57467196456786, -2.29115216667741, -6.82361882016302,
                         -3.23314126063932, 5.03687573600115, -7.7412062869873, -0.455529780242026, -6.0503439372704};
  std::cout << G.inverse() << '\n';

  return 0;
}