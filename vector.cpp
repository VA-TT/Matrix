#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <initializer_list>
#include <cmath>
#include <iomanip>

template <typename T>
class Vector
{
private:
    std::vector<T> m_elements;

public:
    // Default Constructor
    Vector() = default;
    
    // Constructor với kích thước n, tất cả phần tử = 0
    explicit Vector(std::size_t n) : m_elements(n, T{}) {}
    
    // Constructor với kích thước n, tất cả phần tử = value
    Vector(std::size_t n, const T& value) : m_elements(n, value) {}
    
    // Constructor với initializer_list
    Vector(std::initializer_list<T> list) : m_elements(list) {}
    
    // Constructor 3D (x, y, z)
    Vector(const T& x, const T& y, const T& z) : m_elements{x, y, z} {}
    
    // Getters
    std::size_t size() const { return m_elements.size(); }
    bool empty() const { return m_elements.empty(); }
    
    // Element access
    T& operator[](std::size_t i) 
    { 
        assert(i < size());
        return m_elements[i]; 
    }
    
    const T& operator[](std::size_t i) const 
    { 
        assert(i < size());
        return m_elements[i]; 
    }
    
    T& at(std::size_t i) { return m_elements.at(i); }
    const T& at(std::size_t i) const { return m_elements.at(i); }
    
    // Iterators
    auto begin() { return m_elements.begin(); }
    auto end() { return m_elements.end(); }
    auto begin() const { return m_elements.begin(); }
    auto end() const { return m_elements.end(); }
    
    // Resize vector
    void resize(std::size_t n) { m_elements.resize(n); }
    void resize(std::size_t n, const T& value) { m_elements.resize(n, value); }
    
    // Add element
    void push_back(const T& value) { m_elements.push_back(value); }
    
    // Print
    friend std::ostream& operator<<(std::ostream& os, const Vector<T>& v)
    {
        os << "[";
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            os << v[i];
            if (i < v.size() - 1) os << ", ";
        }
        os << "]";
        return os;
    }
};

// Vector addition
template <typename T>
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2)
{   
    if (v1.size() != v2.size())
        throw std::invalid_argument("Vectors must have the same dimension.");
    
    Vector<T> result(v1.size());
    for (std::size_t i = 0; i < v1.size(); ++i)
        result[i] = v1[i] + v2[i];  
    return result;
}

// Scalar multiplication (scalar * vector)
template <typename T>
Vector<T> operator*(const T& k, const Vector<T>& v)
{   
    Vector<T> result(v.size());
    for (std::size_t i = 0; i < v.size(); ++i)
        result[i] = k * v[i];  
    return result;
}

// Scalar multiplication (vector * scalar)
template <typename T>
Vector<T> operator*(const Vector<T>& v, const T& k)
{    
    return k * v;
}

// Unary minus
template <typename T>
Vector<T> operator-(const Vector<T>& v)
{    
    return T{-1} * v;
}

// Vector subtraction
template <typename T>
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2)
{    
    return v1 + (-v2);
}

// Dot product
template <typename T>
T dotProduct(const Vector<T>& v1, const Vector<T>& v2)
{
    if (v1.size() != v2.size())
        throw std::invalid_argument("Vectors must have the same dimension.");

    T result{};
    for (std::size_t i = 0; i < v1.size(); ++i)
        result += v1[i] * v2[i];  
    return result;
}

// Cross product (chỉ cho vector 3D)
template <typename T>
Vector<T> crossProduct(const Vector<T>& v1, const Vector<T>& v2)
{
    if (v1.size() != 3 || v2.size() != 3)
        throw std::invalid_argument("Cross product is only defined for 3D vectors.");
    
    Vector<T> result(3);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

// Vector magnitude/norm
template <typename T>
T magnitude(const Vector<T>& v)
{
    return std::sqrt(dotProduct(v, v));
}

// Unit vector
template <typename T>
Vector<T> normalize(const Vector<T>& v)
{
    T mag = magnitude(v);
    if (mag == T{})
        throw std::invalid_argument("Cannot normalize zero vector.");
    return v * (T{1} / mag);
}