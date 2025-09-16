#include <array>
#include <iostream>
#include <algorithm>
#include <cassert>
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
      result(i, i) = 1;
    return result;
  }
  void resetIdentity()
  {
    (*this) = Matrix::identity();
  }

  // Accesing the elements in the array with one parameter (i)
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

  // Accesing the elements in the matrices with two parameters (i,j)
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

  // Reference to a whole row(i) or column(j)
  std::array<std::reference_wrapper<T>, nCols> referenceRow(std::size_t i)
  {
    assert(i < nRows);
    std::array<std::reference_wrapper<T>, nCols> row_refs{};
    for (std::size_t j{0}; j < nCols; ++j)
    {
      row_refs[j] = std::ref((*this)(i, j));
    }
    return row_refs;
  }

  std::array<std::reference_wrapper<T>, nRows> referenceCol(std::size_t j)
  {
    assert(j < nCols);
    std::array<std::reference_wrapper<T>, nRows> col_refs{};
    for (std::size_t i{0}; i < nRows; ++i)
    {
      col_refs[i] = std::ref((*this)(i, j));
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
    if (m1.getCols() != m2.getCols() || m1.getRows() != m2.getRows())
    {
      return false;
    }
    for (std::size_t i{0}; i < m1.length(); i++)
    {
      if (m1[i] != m2[i])
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

  // Swap 2 rows
  void swapRows(std::size_t i1, std::size_t i2)
  {
    for (std::size_t j{0}; j < nCols; ++j)
    {
      std::swap((*this)(i1, j), (*this)(i2, j);)
    }
  }

  // Find the row with the maximum number of elements that aren't 0
  std::size_t LeastZeroRow(double tolerance = 5e-4)
  {
    std::array<int, nRows> exceed0{};
    std::size_t index = 0;
    for (std::size_t i{0}; i < nRows; i++)
    {
      for (auto &element : (*this).referenceRow(i))
      {
        if (element.get() < tolerance)
          ++exceed0[i];
      }
    }
    for (std::size_t i{0}; i < nRows; i++)
    {
      if (exceed0[i] > exceed0[index])
      {
        index = i;
      }
    }
    // Tim index chua gia tri max trong array exceed0
    return index;
  }

  // Multiply a const to a row

  // Adding a row multiplied with a const to another row

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
    // check det
    Matrix<T, nRows, nCols> A(*this);
    Matrix<T, nRows, nCols> I = Matrix<T, nRows, nCols>::identity();

    for (std::size_t i = 0; i < nRows; ++i)
    {
      // Find pivot
      T pivot = A(i, i);
      std::size_t pivotRow = i;
      for (std::size_t row = i + 1; row < nRows; ++row)
      {
        if (std::abs(A(row, i)) > std::abs(pivot))
        {
          pivot = A(row, i);
          pivotRow = row;
        }
      }
      if (pivot == 0)
      {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
      }
      // Swap rows if needed
      if (pivotRow != i)
      {
        for (std::size_t col = 0; col < nCols; ++col)
        {
          std::swap(A(i, col), A(pivotRow, col));
          std::swap(I(i, col), I(pivotRow, col));
        }
      }
      // Normalize pivot row
      T invPivot = 1 / A(i, i);
      for (std::size_t col = 0; col < nCols; ++col)
      {
        A(i, col) *= invPivot;
        I(i, col) *= invPivot;
      }
      // Eliminate other rows
      for (std::size_t row = 0; row < nRows; ++row)
      {
        if (row == i)
          continue;
        T factor = A(row, i);
        for (std::size_t col = 0; col < nCols; ++col)
        {
          A(row, col) -= factor * A(i, col);
          I(row, col) -= factor * I(i, col);
        }
      }
    }
    return I;
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

template <typename T, std::size_t R1, std::size_t C1, std::size_t R2, std::size_t C2>
Matrix<T, R1, (C1 + C2)> augmentedMatrix(const Matrix<T, R1, C1> &A, const Matrix<T, R2, C2> &B)
{
  assert(R1 == R2 && "The number of rows of both matrices must match to create an augmented matrix.\n");
  Matrix<T, R1, (C1 + C2)> augmentedMatrix{};
  for (std::size_t i = 0; i < R1; ++i)
    for (std::size_t j = 0; j < (C1 + C2); ++j)
    {
      if (j < C1)
      {
        assert(j < A.getCols() && "Error index exceeding matrix's size.\n");
        augmentedMatrix(i, j) = A(i, j);
      }
      else
      {
        assert((j - C1) < B.getCols() && "Error index exceeding matrix's size.\n");
        augmentedMatrix(i, j) = B(i, j - C1);
      }
    }
  return augmentedMatrix;
}

// Separate a matrix to two matrices by a column

// Gaussian elimination to solve Ax = B
template <typename T, std::size_t R1, std::size_t C1>
Matrix<T, R1, 1> gaussianElimination(const Matrix<T, R1, C1> &A, const Matrix<T, R1, 1> &B)
{
  // assert(R1 == R2 && "Not suitable for solving Ax = B.");
  std::size_t augCols{C1 + 1};
  Matrix<T, R1, (augCols)> augmentedMatrix{};
  for (std::size_t i = 0; i < R1; ++i)
    for (std::size_t j = 0; j < C1; ++j)
    {
      augmentedMatrix(i, j) = A(i, j);
      augmentedMatrix(i, C1) = B(i, 0);
    }

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

int main()
{
  // Matrix<int, 2, 3> A{2, 1, -1,
  //                     1, -1, 1};
  Matrix<int, 2, 3> B{4, -2, 1,
                      2, -4, -2};
  auto row_refs = B.referenceRow(1);
  for (auto &ref : row_refs)
  {
    ref.get() = 42; // thay đổi giá trị từng phần tử
  }
  // Matrix<int, 2, 2> C{1, 2,
  //                     2, 1};
  // Matrix<int, 2, 2> D{3, 4,
  //                     4, 3};
  // Matrix<int, 1, 2> E{1, 2};
  // Matrix<int, 2, 2> G{1, 2,
  //                     2, 4};
  // Matrix<int, 2, 2> H{2, 1,
  //                     1, 3};
  // Matrix<int, 2, 2> I{4, 3,
  //                     0, 2};

  // std::cout << B - 2 * A << '\n';
  // std::cout << C + E << '\n';
  // // std::cout << A * C << '\n';
  // std::cout << C * D << '\n';
  // std::cout << C * B << '\n';
  // std::cout << std::boolalpha << (G * H == G * I) << '\n';
  // std::cout << std::boolalpha << (H != I) << '\n';

  // Matrix<int, 3, 3> A{1, 1, 1,
  //                     1, 2, 3,
  //                     1, 3, 4};
  // Matrix<int, 3, 3> D{2, 0, 0,
  //                     0, 3, 0,
  //                     0, 0, 4};
  // std::cout << B.transpose() << '\n';
  // std::cout << (A * D).transpose() << '\n';
  // std::cout << (D.transpose() * A.transpose()) << '\n';

  // Matrix<int, 2, 2> matrix1{1, 1, 1, 0};
  // Matrix<int, 2, 2> matrix2{5, 3, 3, 2};
  // std::cout << matrix1 * matrix2 << '\n';
  // Matrix<int, 3, 3> matrix3{Matrix<int, 3, 3>::identity()};
  // std::cout << std::boolalpha << matrix3.isOrthogonal() << '\n';
  Matrix<int, 3, 3> A{3, -7, -2, -3, 5, 1, 6, -4, 0};
  return 0;
}