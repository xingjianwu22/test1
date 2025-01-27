#pragma once
#include "ray.hpp"
#include "Vector3f.h"
#include <limits>
#include <array>

class Bounds3
{
  public:
    Vector3f pMin, pMax; // two points to specify the bounding box
    Bounds3()
    {
        double minNum = std::numeric_limits<double>::lowest();
        double maxNum = std::numeric_limits<double>::max();
        pMax = Vector3f(minNum, minNum, minNum);
        pMin = Vector3f(maxNum, maxNum, maxNum);
    }
    Bounds3(const Vector3f p) : pMin(p), pMax(p) {}
    Bounds3(const Vector3f p1, const Vector3f p2)
    {
        pMin = Vector3f(fmin(p1.x(), p2.x()), fmin(p1.y(), p2.y()), fmin(p1.z(), p2.z()));
        pMax = Vector3f(fmax(p1.x(), p2.x()), fmax(p1.y(), p2.y()), fmax(p1.z(), p2.z()));
    }

    Vector3f Diagonal() const { return pMax - pMin; }
    int maxExtent() const
    {
        Vector3f d = Diagonal();
        if (d.x() > d.y() && d.x() > d.z())
            return 0;
        else if (d.y() > d.z())
            return 1;
        else
            return 2;
    }

    double SurfaceArea() const
    {
        Vector3f d = Diagonal();
        return 2 * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
    }

    Vector3f Centroid() { return 0.5 * pMin + 0.5 * pMax; }
    Bounds3 Intersect(const Bounds3& b)
    {
        return Bounds3(Vector3f(fmax(pMin.x(), b.pMin.x()), fmax(pMin.y(), b.pMin.y()),
                                fmax(pMin.z(), b.pMin.z())),
                       Vector3f(fmin(pMax.x(), b.pMax.x()), fmin(pMax.y(), b.pMax.y()),
                                fmin(pMax.z(), b.pMax.z())));
    }

    Vector3f Offset(const Vector3f& p) const
    {
        Vector3f o = p - pMin;
        if (pMax.x() > pMin.x())
            o.x() /= pMax.x() - pMin.x();
        if (pMax.y() > pMin.y())
            o.y() /= pMax.y() - pMin.y();
        if (pMax.z() > pMin.z())
            o.z() /= pMax.z() - pMin.z();
        return o;
    }

    bool Overlaps(const Bounds3& b1, const Bounds3& b2)
    {
        bool x = (b1.pMax.x() >= b2.pMin.x()) && (b1.pMin.x() <= b2.pMax.x());
        bool y = (b1.pMax.y() >= b2.pMin.y()) && (b1.pMin.y() <= b2.pMax.y());
        bool z = (b1.pMax.z() >= b2.pMin.z()) && (b1.pMin.z() <= b2.pMax.z());
        return (x && y && z);
    }

    bool Inside(const Vector3f& p, const Bounds3& b)
    {
        return (p.x() >= b.pMin.x() && p.x() <= b.pMax.x() && p.y() >= b.pMin.y() &&
                p.y() <= b.pMax.y() && p.z() >= b.pMin.z() && p.z() <= b.pMax.z());
    }
    inline const Vector3f& operator[](int i) const
    {
        return (i == 0) ? pMin : pMax;
    }

    inline bool IntersectP(const Ray& ray, const Vector3f& invDir,
                           const std::array<int, 3>& dirisNeg) const;
};



inline bool Bounds3::IntersectP(const Ray& ray, const Vector3f& invDir,
                                const std::array<int, 3>& dirIsNeg) const
{
    // ray.orig.x + dir.x * tx_min = pmin.x
    // tx_min = (pmin.x - orig.x) * invdir.x
    double tx_min = (pMin.m_elements[0] - ray.origin.m_elements[0]) * invDir.m_elements[0];
    double tx_max = (pMax.m_elements[0] - ray.origin.m_elements[0]) * invDir.m_elements[0];
    if (!dirIsNeg[0])
        std::swap(tx_min, tx_max);

    double ty_min = (pMin.m_elements[1] - ray.origin.m_elements[1]) * invDir.m_elements[1];
    double ty_max = (pMax.m_elements[1] - ray.origin.m_elements[1]) * invDir.m_elements[1];
    if (!dirIsNeg[1])
        std::swap(ty_min, ty_max);
    
    double tz_min = (pMin.m_elements[2] - ray.origin.m_elements[2]) * invDir.m_elements[2];
    double tz_max = (pMax.m_elements[2] - ray.origin.m_elements[2]) * invDir.m_elements[2];
    if (!dirIsNeg[2])
        std::swap(tz_min, tz_max);
    // 更新 tenter，取最大值，因为我们需要射线在三个维度上同时进入包围框的最晚时间。
    // 更新 texit，取最小值，因为我们需要射线在三个维度上存在离开包围框的最早时间
    double t_enter = std::max(tx_min, std::max(ty_min, tz_min));
    double t_exit = std::min(tx_max, std::min(ty_max, tz_max));
    // tenter <= texit: 射线进入时间必须早于射线退出时间，这表示射线确实穿过了包围框。
    // texit >= 0: 射线退出时间必须大于等于零，这表示包围框在射线的正向路径上。
    return t_enter <= t_exit && t_exit >= 0;
}

inline Bounds3 Union(const Bounds3& b1, const Bounds3& b2)
{
    Bounds3 ret;
    ret.pMin = Vector3f::Min(b1.pMin, b2.pMin);
    ret.pMax = Vector3f::Max(b1.pMax, b2.pMax);
    return ret;
}

inline Bounds3 Union(const Bounds3& b, const Vector3f& p)
{
    Bounds3 ret;
    ret.pMin = Vector3f::Min(b.pMin, p);
    ret.pMax = Vector3f::Max(b.pMax, p);
    return ret;
}
