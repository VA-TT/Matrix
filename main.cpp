#include <array>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iomanip>

template <typename T, std::size_t nRows, std::size_t nCols>
class Matrix
{
private:
  std::array<T, nRows * nCols> m_elements{};

public:
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
  T &operator()(std::size_t i, std::size_t j)
  {
    assert(i < nRows && j < nCols);
    return m_elements[i * nCols + j];
  }
  const T &operator()(std::size_t i, std::size_t j) const
  {
    assert(i < nRows && j < nCols);
    return m_elements[i * nCols + j];
  }

  constexpr std::size_t cols() const { return nCols; }
  constexpr std::size_t rows() const { return nRows; }
  constexpr int length() const { return nRows * nCols; }

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

  // Matrix algebra
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
    assert(m1.cols() == m2.cols() && m1.rows() == m2.rows() && "Unable to perform matrix addition/substraction.");
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
    if (m1.cols() != m2.cols() || m1.rows() != m2.rows())
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

  // Matrix<T, nCols, nRows> inverse() const {}

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
};

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

template <typename T, std::size_t R1, std::size_t C1, std::size_t R2, std::size_t C2>
bool arePairOrthogonal(const Matrix<T, R1, C1> &m1, const Matrix<T, R2, C2> &m2)
{
  if ((m1 != Matrix<T, R1, C1>::zero()) && (m2 != Matrix<T, R2, C2>::zero()))
  {
    return ((m1.transpose() * m2) == Matrix<T, R1, C2>::identity());
  }
  return false; // Both matrices must not be 0
}

template <typename T, std::size_t nRows, std::size_t nCols>
T trace(const Matrix<T, nRows, nCols> &m)
{
  static_assert(nRows == nCols, "Trace chỉ áp dụng cho ma trận vuông");
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
  Matrix<int, 3, 3> B{-7, 5, 2}

  return 0;
}