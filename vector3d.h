#ifndef Vector3d_H
#define Vector3d_H

// class : vector double [3]
class Vector3d
{
    public:
        double x, y, z;

        // Constructor
        Vector3d();
        Vector3d(double x, double y, double z);

        // Assignment operator
        Vector3d& operator=(const Vector3d& v);

        // Unary operator
        Vector3d& operator+=(const Vector3d& v);
        Vector3d& operator-=(const Vector3d& v);
        Vector3d& operator*=(double k);
        Vector3d& operator/=(double k);
        Vector3d operator+() const;
        Vector3d operator-() const;

        // Subscript operator
        double& operator[](int i);

        // Comparioson operator
        bool operator==(const Vector3d& v) const;
        bool operator!=(const Vector3d& v) const;

        // Normalize
        double norm() const;
        void normalize();

        void initialize();
        void set_value(double xx, double yy, double zz);
};

/* - vector3d + vector3d - */
Vector3d operator+(const Vector3d& u, const Vector3d& v);

/* - vector3d - vector3d - */
Vector3d operator-(const Vector3d& u, const Vector3d& v);

/* - double   * vector3d - */
Vector3d operator*(double k, const Vector3d& v);

/* - vector3d * double   - */
Vector3d operator*(const Vector3d& v, double k);

/* - vector3d / double   - */
Vector3d operator/(const Vector3d& v, double k);

/* - inner product vector3d * vector3d - */
double operator*(const Vector3d& u, const Vector3d& v);

/* - outer product vector3d % vector3d - */
Vector3d operator%(const Vector3d& u, const Vector3d& v);

#include <cmath>
#include <iostream>
using namespace std;

inline Vector3d::Vector3d() { x = y = z = 0.0; }
inline Vector3d::Vector3d(double x, double y, double z)
{
    this->x = x;  this->y = y;  this->z = z;
}

inline Vector3d& Vector3d::operator=(const Vector3d& v)
{
    if (this != &v) {
        this->x = v.x;  this->y = v.y;  this->z = v.z;
    }
    return *this;
}

inline Vector3d& Vector3d::operator+=(const Vector3d& v)
{
    this->x += v.x;  this->y += v.y;  this->z += v.z;
    return *this;
}

inline Vector3d& Vector3d::operator-=(const Vector3d& v)
{
    this->x -= v.x;  this->y -= v.y;  this->z -= v.z;
    return *this;
}

inline Vector3d& Vector3d::operator*=(double k)
{
    this->x *= k;  this->y *= k;  this->z *= k;
    return *this;
}

inline Vector3d& Vector3d::operator/=(double k)
{
    this->x /= k;  this->y /= k;  this->z /= k;
    return *this;
}

inline Vector3d Vector3d::operator+()const { return *this; }
inline Vector3d Vector3d::operator-()const { return Vector3d(-x, -y, -z); }

inline double& Vector3d::operator[](int i)
{
    if (i == 0) {
        return x;
    } else if (i == 1) {
        return y;
    } else if (i == 2) {
        return z;
    } else {
        return x;
    }
}

inline bool Vector3d::operator==(const Vector3d& v) const
{
    return (x == v.x)&&(y == v.y)&&(z == v.z);
}

inline bool Vector3d::operator!=(const Vector3d& v) const
{
    return !(*this == v);
}

inline double Vector3d::norm() const
{
    return sqrt(x * x + y * y + z * z);
}

inline void Vector3d::normalize()
{
    *this /= norm();
}

inline void Vector3d::initialize() { x = y = z = 0.0;}

inline void Vector3d::set_value(double xx, double yy, double zz)
{
    this->x = xx;   this->y = yy;   this->z = zz;
}

/* - vector3d + vector3d - */
inline Vector3d operator+(const Vector3d& u, const Vector3d& v)
{
    Vector3d w;
    w.x = u.x + v.x;  w.y = u.y + v.y;  w.z = u.z + v.z;
    return w;
}

/* - vector3d - vector3d - */
inline Vector3d operator-(const Vector3d& u, const Vector3d& v)
{
    Vector3d w;
    w.x = u.x - v.x;  w.y = u.y - v.y;  w.z = u.z - v.z;
    return w;
}

#endif
