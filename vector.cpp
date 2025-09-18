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

template <typename T>
class Vector
{
private:
    T m_x{};
    T m_y{};
    T m_z{};
    std::array<T, 3> m_vector{m_x, m_y, m_z};

public:
} Vector operator+(const Vector &v1, const Vector &v2)
{
    Vector v{};
    for (T element : m_vector)
    {
        element =
    }
}