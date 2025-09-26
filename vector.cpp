#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <initializer_list>
#include <cmath>
#include <iomanip>
#include "comparison.hpp"   //Approximative Comparsion

using Index = std::ptrdiff_t;  // typedef

template <typename T>
class Vect3d
{
private:
    std::vector<T> m_elements{};

public:
    // Default Constructors
    Vect3d() = default;
    Vect3d(const Vect3d &) = default;
    Vect3d(Vect3d &&) = default;
    Vect3d &operator=(const Vect3d &) = default;  // Sửa: constVect3dx → const Vect3d
    Vect3d &operator=(Vect3d &&) = default;
    ~Vect3d() = default;
    
    // Constructor with length n of elements which are 0
    explicit Vect3d(std::size_t n) : m_elements(n, T{}) {}
    
    // // Constructor with length n of elements which are value
    // Vect3d(std::size_t n, const T& value) : m_elements(n, value) {}
    
    // Constructor với initializer_list
    Vect3d(std::initializer_list<T> list) : m_elements(list) {}
    
    // // Constructor 3D (x, y, z)
    // Vect3d(const T& x, const T& y, const T& z) : m_elements{x, y, z} {}
    
    // Signed indexing với assert
    auto& operator[](Index i) { 
        assert(i >= 0 && "Negative index not allowed");
        assert(static_cast<std::size_t>(i) < m_elements.size() && "Index out of bounds");
        return m_elements.data()[static_cast<std::size_t>(i)]; 
    }
    
    const auto& operator[](Index i) const { 
        assert(i >= 0 && "Negative index not allowed");
        assert(static_cast<std::size_t>(i) < m_elements.size() && "Index out of bounds");
        return m_elements.data()[static_cast<std::size_t>(i)]; 
    }
    

    
    // Size functions
    // constexpr std::size_t size() const noexcept { return m_elements.size(); }
    constexpr Index size() const noexcept { 
        return static_cast<Index>(m_elements.size()); 
    }
    
    
    T& at(Index i) { 
        if (i < 0) throw std::out_of_range("Negative index not allowed");
        return m_elements.at(static_cast<std::size_t>(i)); 
    }
    
    const T& at(Index i) const { 
        if (i < 0) throw std::out_of_range("Negative index not allowed");
        return m_elements.at(static_cast<std::size_t>(i)); 
    }
    
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
    friend std::ostream& operator<<(std::ostream& os, const Vect3d<T>& v)
    {
        os << "[";
        for (Index i = 0; i < v.size(); ++i)
        {
            os << v[i];
            if (i < v.size() - 1) os << ", ";
        }
        os << "]";
        return os;
    }

    //Projection on another vector
    Vect3d<T> projection (const Vect3d& other) const
    {
        return normalize(other)*dotProduct(*this, normalize(other));
    }
};


// Vect3d operators
template <typename T>
bool operator==(const Vect3d<T>& v1, const Vect3d<T>& v2)
{   
    if (v1.size() != v2.size())
        return false;
    for (Index i = 0; i < v1.size(); ++i)
    { 
        if (!approximatelyEqualAbsRel(v1[i], v2[i]))
            return false;
    }
    return true;
}

template <typename T>
bool operator!=(const Vect3d<T>& v1, const Vect3d<T>& v2)
{   
    return !(v1 == v2);
}

template <typename T>
Vect3d<T> operator+(const Vect3d<T>& v1, const Vect3d<T>& v2)
{   
    if (v1.size() != v2.size())
        throw std::invalid_argument("Vectors must have the same dimension.");
    
    Vect3d<T> result(v1.size());
    for (Index i = 0; i < v1.size(); ++i)
        result[i] = v1[i] + v2[i];  
    return result;
}

// Scalar multiplication (scalar * vector)
template <typename T>
Vect3d<T> operator*(const T& k, const Vect3d<T>& v)
{   
    Vect3d<T> result(v.size());
    for (Index i = 0; i < v.size(); ++i)
        result[i] = k * v[i];  
    return result;
}

// Scalar multiplication (vector * scalar)
template <typename T>
Vect3d<T> operator*(const Vect3d<T>& v, const T& k)
{    
    return k * v;
}

// Unary minus
template <typename T>
Vect3d<T> operator-(const Vect3d<T>& v)
{    
    return T{-1} * v;
}

// Vect3d subtraction
template <typename T>
Vect3d<T> operator-(const Vect3d<T>& v1, const Vect3d<T>& v2)
{    
    return v1 + (-v2);
}

// Dot product
template <typename T>
T dotProduct(const Vect3d<T>& v1, const Vect3d<T>& v2)
{
    if (v1.size() != v2.size())
        throw std::invalid_argument("Vectors must have the same dimension.");

    T result{};
    for (Index i = 0; i < v1.size(); ++i)
        result += v1[i] * v2[i];  
    return result;
}

// Cross product (chỉ cho vector 3D)
template <typename T>
Vect3d<T> crossProduct(const Vect3d<T>& v1, const Vect3d<T>& v2)
{
    if (v1.size() != 3 || v2.size() != 3)
        throw std::invalid_argument("Cross product is only defined for 3D vectors.");
    
    Vect3d<T> result(3);
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

// Vect3d magnitude/norm
template <typename T>
T magnitude(const Vect3d<T>& v)
{
    return std::sqrt(dotProduct(v, v));
}

// Unit vector
template <typename T>
Vect3d<T> normalize(const Vect3d<T>& v)
{
    T mag = magnitude(v);
    if (mag == T{})
        throw std::invalid_argument("Cannot normalize zero vector.");
    return v * (T{1} / mag);
}

// Unit vector
template <typename T>
T angleRad(const Vect3d<T>& v1, const Vect3d<T>& v2)
{
    return std::acos(dotProduct(v1,v2) / (magnitude(v1)*magnitude(v2)));
}

template <typename T>
T angleDegree(const Vect3d<T>& v1, const Vect3d<T>& v2)
{
    const T PI = T{3.14159265358979323846};
    return angleRad(v1, v2) * T{180} / PI;
}

template <typename T>
bool isPerpendicular(const Vect3d<T>& v1, const Vect3d<T>& v2)
{
    return approximatelyEqualAbsRel(dotProduct(v1,v2), T{});
}

template <typename T>
bool isParallel(const Vect3d<T>& v1, const Vect3d<T>& v2)
{
    return approximatelyEqualAbsRel(magnitude(crossProduct(v1,v2)), T{});
}


int main() {
    std::cout << "=== VECTOR CLASS TEST SUITE ===" << std::endl;
    std::cout << std::endl;
    
    // Test 1: Constructors
    std::cout << "1. Testing Constructors:" << std::endl;
    Vect3d<int> v1;                    // Default constructor
    Vect3d<int> v2(5);                 // Size constructor (5 zeros)
    Vect3d<int> v3{1, 2, 3, 4, 5};     // Initializer list
    
    std::cout << "v1 (default): " << v1 << " (size: " << v1.size() << ")" << std::endl;
    std::cout << "v2(5): " << v2 << " (size: " << v2.size() << ")" << std::endl;
    std::cout << "v3{1,2,3,4,5}: " << v3 << " (size: " << v3.size() << ")" << std::endl;
    std::cout << std::endl;
    
    // Test 2: Indexing
    std::cout << "2. Testing Indexing:" << std::endl;
    Vect3d<double> v4{10.5, 20.3, 30.7, 40.1};
    std::cout << "v4: " << v4 << std::endl;
    std::cout << "v4[0] = " << v4[0] << std::endl;
    std::cout << "v4[2] = " << v4[2] << std::endl;
    
    // Modify element
    v4[1] = 99.9;
    std::cout << "After v4[1] = 99.9: " << v4 << std::endl;
    
    // Test signed indexing
    Index idx = 3;
    std::cout << "Using Index idx=3: v4[idx] = " << v4[idx] << std::endl;
    std::cout << std::endl;
    
    // Test 3: Vect3d Operations
    std::cout << "3. Testing Vect3d Operations:" << std::endl;
    Vect3d<double> va{1.0, 2.0, 3.0};
    Vect3d<double> vb{4.0, 5.0, 6.0};
    
    std::cout << "va = " << va << std::endl;
    std::cout << "vb = " << vb << std::endl;
    
    auto vc = va + vb;
    std::cout << "va + vb = " << vc << std::endl;
    
    auto vd = vb - va;
    std::cout << "vb - va = " << vd << std::endl;
    
    auto ve = 2.5 * va;
    std::cout << "2.5 * va = " << ve << std::endl;
    
    auto vf = va * 3.0;
    std::cout << "va * 3.0 = " << vf << std::endl;
    
    auto vg = -va;
    std::cout << "-va = " << vg << std::endl;
    std::cout << std::endl;
    
    // Test 4: Vect3d Math
    std::cout << "4. Testing Vect3d Math:" << std::endl;
    Vect3d<double> u1{3.0, 4.0, 0.0};
    Vect3d<double> u2{1.0, 0.0, 0.0};
    
    std::cout << "u1 = " << u1 << std::endl;
    std::cout << "u2 = " << u2 << std::endl;
    
    double dot = dotProduct(u1, u2);
    std::cout << "Dot product u1·u2 = " << dot << std::endl;
    
    auto cross = crossProduct(u1, u2);
    std::cout << "Cross product u1×u2 = " << cross << std::endl;
    
    double mag1 = magnitude(u1);
    double mag2 = magnitude(u2);
    std::cout << "Magnitude |u1| = " << mag1 << std::endl;
    std::cout << "Magnitude |u2| = " << mag2 << std::endl;
    
    auto unit1 = normalize(u1);
    auto unit2 = normalize(u2);
    std::cout << "Unit vector of u1 = " << unit1 << std::endl;
    std::cout << "Unit vector of u2 = " << unit2 << std::endl;
    std::cout << "Magnitude of unit u1 = " << magnitude(unit1) << std::endl;
    std::cout << std::endl;
    
    // Test 5: Resize and Push Back
    std::cout << "5. Testing Resize and Push Back:" << std::endl;
    Vect3d<int> vr{10, 20, 30};
    std::cout << "Original vr: " << vr << " (size: " << vr.size() << ")" << std::endl;
    
    vr.resize(6);
    std::cout << "After resize(6): " << vr << " (size: " << vr.size() << ")" << std::endl;
    
    vr.resize(8, 77);
    std::cout << "After resize(8, 77): " << vr << " (size: " << vr.size() << ")" << std::endl;
    
    vr.push_back(100);
    std::cout << "After push_back(100): " << vr << " (size: " << vr.size() << ")" << std::endl;
    std::cout << std::endl;
    
    // Test 6: Edge Cases and Error Handling
    std::cout << "6. Testing Edge Cases:" << std::endl;
    
    try {
        Vect3d<int> v_small{1, 2};
        Vect3d<int> v_big{1, 2, 3, 4};
        std::cout << "Trying to add vectors of different sizes..." << std::endl;
        auto result = v_small + v_big;
    } catch (const std::exception& e) {
        std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
    }
    
    try {
        Vect3d<double> zero_vec{0.0, 0.0, 0.0};
        std::cout << "Trying to normalize zero vector..." << std::endl;
        auto unit = normalize(zero_vec);
    } catch (const std::exception& e) {
        std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
    }
    
    try {
        Vect3d<int> v2d_1{1, 2};
        Vect3d<int> v2d_2{3, 4};
        std::cout << "Trying cross product on 2D vectors..." << std::endl;
        auto cross_2d = crossProduct(v2d_1, v2d_2);
    } catch (const std::exception& e) {
        std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
    }
    
    std::cout << std::endl;
    std::cout << "=== ALL TESTS COMPLETED SUCCESSFULLY! ===" << std::endl;
    
    // Test 7: Projection
    std::cout << "7. Testing Vector Projection:" << std::endl;

    Vect3d<double> proj_a{3.0, 4.0, 0.0};      // Vector a
    Vect3d<double> proj_b{1.0, 0.0, 0.0};      // Unit vector along x-axis
    Vect3d<double> proj_c{2.0, 2.0, 0.0};      // 45-degree vector

    std::cout << "Vector proj_a = " << proj_a << std::endl;
    std::cout << "Vector proj_b = " << proj_b << std::endl;
    std::cout << "Vector proj_c = " << proj_c << std::endl;

    // Vector projection
    auto vector_proj_ab = proj_a.projection(proj_b);
    auto vector_proj_ac = proj_a.projection(proj_c);
    std::cout << "Projection of proj_a onto proj_b: " << vector_proj_ab << std::endl;
    std::cout << "Projection of proj_a onto proj_c: " << vector_proj_ac << std::endl;

    std::cout << std::endl;

    // Test 8: Angle Calculation
    std::cout << "8. Testing Angle Calculation:" << std::endl;

    Vect3d<double> angle_v1{1.0, 0.0, 0.0};  // Unit vector along x-axis
    Vect3d<double> angle_v2{0.0, 1.0, 0.0};  // Unit vector along y-axis  
    Vect3d<double> angle_v3{1.0, 1.0, 0.0};  // 45-degree vector

    std::cout << "angle_v1 = " << angle_v1 << std::endl;
    std::cout << "angle_v2 = " << angle_v2 << std::endl;
    std::cout << "angle_v3 = " << angle_v3 << std::endl;

    // Test 90 degree angle
    double angle12_rad = angleRad(angle_v1, angle_v2);
    double angle12_deg = angleDegree(angle_v1, angle_v2);
    std::cout << "Angle between angle_v1 and angle_v2: " << angle12_rad << " rad = " << angle12_deg << " deg" << std::endl;

    // Test 45 degree angle
    double angle13_rad = angleRad(angle_v1, angle_v3);
    double angle13_deg = angleDegree(angle_v1, angle_v3);
    std::cout << "Angle between angle_v1 and angle_v3: " << angle13_rad << " rad = " << angle13_deg << " deg" << std::endl;

    std::cout << std::endl;

    // Test 9: Perpendicular and Parallel
    std::cout << "9. Testing Perpendicular and Parallel:" << std::endl;

    Vect3d<double> perp_x{1.0, 0.0, 0.0};      // Unit vector along x
    Vect3d<double> perp_y{0.0, 1.0, 0.0};      // Unit vector along y  
    Vect3d<double> para_x{2.0, 0.0, 0.0};      // Parallel to perp_x

    std::cout << "perp_x = " << perp_x << std::endl;
    std::cout << "perp_y = " << perp_y << std::endl;
    std::cout << "para_x = " << para_x << std::endl;

    // Test perpendicular (should be true)
    bool is_perp = isPerpendicular(perp_x, perp_y);
    std::cout << "isPerpendicular(perp_x, perp_y): " << is_perp << std::endl;

    // Test parallel (should be true)
    bool is_para = isParallel(perp_x, para_x);
    std::cout << "isParallel(perp_x, para_x): " << is_para << std::endl;

    // Show verification values
    std::cout << "Dot product perp_x·perp_y = " << dotProduct(perp_x, perp_y) << std::endl;
    std::cout << "Cross product magnitude |perp_x×para_x| = " << magnitude(crossProduct(perp_x, para_x)) << std::endl;

    std::cout << std::endl;

    return 0;
}