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
    std::vector<T, 3> m_vector{m_x, m_y, m_z};

public:
} ;
template <typename T>
Vector operator+(const Vector &v1, const Vector &v2)
{   
    static_assert(v1.size() == v2.size(), "Vectors don't have the same dimension.");
    Vector result{};
    for (std::size_t i{0}; i < v1.size(); ++v)
        result[i] = v1[i]+v2[i];  
    return result;
}

template <typename T>
Vector operator*( T k, const Vector &v)
{   
     Vector result{};
    for (std::size_t i{0}; i < v.size(); ++v)
        result[i] = k * v[i];  
    return result;
}

template <typename T>
Vector operator*(const Vector &v, T k)
{    
    return k * v;
}

template <typename T>
Vector operator-(const Vector &v)
{    
    return -1 * v;
}

template <typename T>
Vector operator-(const Vector &v1, const Vector& v2)
{    
    return v1 + (-v2);
}

Vector dotProduct(){}
Vector crossProduct(){}