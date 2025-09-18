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
  std::array<T, nCols> row(std::size_t i) const
  {
    assert(i < nRows);
    std::array<T, nCols> row_i{};
    for (std::size_t j{0}; j < nCols; ++j)
      row_i[j] = (*this)(i, j);
    return row_i;
  }

  std::array<T, nRows> col(std::size_t j) const
  {
    assert(j < nCols);
    std::array<T, nRows> col_j{};
    for (std::size_t i{0}; i < nRows; ++i)
      col_j[i] = (*this)(i, j);
    return col_j;
  }

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
    for (std::size_t i{0}; i < m1.length(); i++)
    {
      if (std::abs(m1[i] - m2[i]) > tolerance)
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

  // Separate a matrix to two matrices by a column
  void splitByColumn(std::size_t colPosition, Matrix &A, Matrix &B)
  {
    assert(A.getRows() == nRows && B.getRows() == nRows && A.getCols() == colPosition && B.getCols() == (nCols - colPosition) && "sum of 2 matrices's size must match the augmented matrix.\n");
    for (std::size_t i = 0; i < nRows; ++i)
      for (std::size_t j = 0; j < nCols; ++j)
      {
        if (j < colPosition)
        {
          assert(j < A.getCols() && "Error index exceeding matrix's size.\n");
          A(i, j) = (*this)(i, j);
        }
        else
        {
          assert(j - colPosition < B.getCols() && "Error index exceeding matrix's size.\n");
          B(i, j - colPosition) = (*this)(i, j);
        }
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

  // Find the row with the max element at a given column
  int indexRowMax(std::size_t j)
  {
    std::array<T, nRows> col_j = (*this).col(j);
    T element_max{};
    std::size_t index{0};
    for (std::size_t i{0}; i < nRows; i++)
    {
      if (col_j[i] > element_max)
      {
        element_max = col_j[i];
        index = i;
      }
    }
    // Tim index chua gia tri max trong array exceed0
    return static_cast<int>(index);
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

  // Inverse Matrix: To be implemented
  Matrix<T, nCols, nRows> inverse() const
  {
    static_assert(nRows == nCols, "Inverse only defined for square matrices.");
  }
  // Matrix<T, nCols, nRows> inverse() const
  // {
  //   static_assert(nRows == nCols, "Inverse only defined for square matrices.");
  //   // check det
  //   Matrix<T, nRows, nCols> A(*this);
  //   Matrix<T, nRows, nCols> I = Matrix<T, nRows, nCols>::identity();

  //   for (std::size_t i = 0; i < nRows; ++i)
  //   {
  //     // Find pivot
  //     T pivot = A(i, i);
  //     std::size_t pivotRow = i;
  //     for (std::size_t row = i + 1; row < nRows; ++row)
  //     {
  //       if (std::abs(A(row, i)) > std::abs(pivot))
  //       {
  //         pivot = A(row, i);
  //         pivotRow = row;
  //       }
  //     }
  //     if (pivot == 0)
  //     {
  //       throw std::runtime_error("Matrix is singular and cannot be inverted.");
  //     }
  //     // Swap rows if needed
  //     if (pivotRow != i)
  //     {
  //       for (std::size_t col = 0; col < nCols; ++col)
  //       {
  //         std::swap(A(i, col), A(pivotRow, col));
  //         std::swap(I(i, col), I(pivotRow, col));
  //       }
  //     }
  //     // Normalize pivot row
  //     T invPivot = 1 / A(i, i);
  //     for (std::size_t col = 0; col < nCols; ++col)
  //     {
  //       A(i, col) *= invPivot;
  //       I(i, col) *= invPivot;
  //     }
  //     // Eliminate other rows
  //     for (std::size_t row = 0; row < nRows; ++row)
  //     {
  //       if (row == i)
  //         continue;
  //       T factor = A(row, i);
  //       for (std::size_t col = 0; col < nCols; ++col)
  //       {
  //         A(row, col) -= factor * A(i, col);
  //         I(row, col) -= factor * I(i, col);
  //       }
  //     }
  //   }
  //   return I;
  // }

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
template <typename T, std::size_t nCols>
std::array<T, nCols> operator*(const T &k, const std::array<T, nCols> &row)
{
  std::array<T, nCols> result{};
  for (std::size_t i = 0; i < nCols; ++i)
    result[i] = row[i] * k;
  return result;
}

template <typename T, std::size_t nCols>
std::array<T, nCols> operator*(const std::array<T, nCols> &row, const T &k)
{
  return k * row;
}
// // Change the value of row directly
// template <typename T>
// std::span<T> operator*(std::span<T> row, const T &k)
// {
//   for (auto &element : row)
//     element *= k;
//   return row;
// }

// template <typename T>
// std::span<T> operator*(const T &k, std::span<T> row)
// {
//   return row * k;
// }

// Adding a row multiplied with a const to another row
template <typename T, std::size_t nCols>
const std::array<T, nCols> operator+(const std::array<T, nCols> &row1, const std::array<T, nCols> &row2)
{
  std::array<T, nCols> result{};
  for (std::size_t j{0}; j < nCols; j++)
  {
    result[j] = row1[j] + row2[j];
  }
  return result;
}

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
  Matrix<int, 2, 3> B{4, -2, 1,
                      2, -4, -2};

  // Test: Thay đổi giá trị từng phần tử của hàng 1
  auto row_refs = B.row(1);
  for (auto &ref : row_refs)
  {
    ref = 42;
  }
  std::cout << "Sau khi thay đổi hàng 1:\n"
            << B << std::endl;

  // Test: Thay đổi giá trị từng phần tử của cột 2
  auto col_refs = B.col(2);
  for (auto &ref : col_refs)
  {
    ref.get() = 99;
  }
  std::cout << "Sau khi thay đổi cột 2:\n"
            << B << std::endl;

  // Test: Cộng, trừ, nhân ma trận
  Matrix<int, 2, 3> A{1, 2, 3,
                      4, 5, 6};
  std::cout << "A + B:\n"
            << A + B << std::endl;
  std::cout << "A - B:\n"
            << A - B << std::endl;
  std::cout << "2 * A:\n"
            << 2 * A << std::endl;

  // Test: Ma trận chuyển vị
  std::cout << "A chuyển vị:\n"
            << A.transpose() << std::endl;

  // Test: So sánh ma trận
  std::cout << std::boolalpha << "A == B? " << (A == B) << std::endl;
  std::cout << std::boolalpha << "A != B? " << (A != B) << std::endl;

  // Test: Trace cho ma trận vuông
  Matrix<int, 3, 3> C{1, 2, 3,
                      4, 5, 6,
                      7, 8, 9};
  std::cout << "Trace(C): " << trace(C) << std::endl;

  return 0;
}