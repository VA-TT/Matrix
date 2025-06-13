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

  // reset zero matrice
  void resetZero()
  {
    for (auto &e : m_elements)
      e = 0;
  }

  void resetIdentity()
  {
    for (std::size_t i = 0; i < this->rows(); ++i)
      for (std::size_t j = 0; j < this->cols(); ++j)
        (*this)(i, j) = (i == j ? 1 : 0);
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

int main()
{
  // Matrix<int, 2, 3> A{2, 1, -1,
  //                     1, -1, 1};
  // Matrix<int, 2, 3> B{4, -2, 1,
  //                     2, -4, -2};
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

  Matrix<int, 3, 3> A{1, 1, 1,
                      1, 2, 3,
                      1, 3, 4};
  Matrix<int, 3, 3> D{2, 0, 0,
                      0, 3, 0,
                      0, 0, 4};
  std::cout << A * D << '\n';
  std::cout << D * A << '\n';
  return 0;
}