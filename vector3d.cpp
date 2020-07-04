#include "vector3d.h"

/* - double   * vector3d - */
Vector3d operator*(double k, const Vector3d& v)
{
    return Vector3d(k * v.x, k * v.y, k * v.z);
}

/* - vector3d * double   - */
Vector3d operator*(const Vector3d& v, double k)
{
    return Vector3d(v.x * k, v.y * k, v.z * k);
}

/* - vector3d / double   - */
Vector3d operator/(const Vector3d& v, double k)
{
    return Vector3d(v.x / k, v.y / k, v.z / k);
}

/* - inner product vector3d * vector3d - */
double operator*(const Vector3d& u, const Vector3d& v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

/* - outer product vector3d % vector3d - */
Vector3d operator%(const Vector3d& u, const Vector3d& v)
{
    Vector3d w;
    w.x = u.y * v.z - u.z * v.y;
    w.y = u.z * v.x - u.x * v.z;
    w.z = u.x * v.y - u.y * v.x;
    return w;
}



