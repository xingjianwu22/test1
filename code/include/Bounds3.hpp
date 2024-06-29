#pragma once
#include "ray.hpp"
#include "Vector3f.h"
#include <array>
#include <cmath>

class Bounds3 
{
public:
    // two points to specify the bounding box
    Vector3f pMin, pMax;

    Bounds3()
    {
        float min_f = std::numeric_limits<float>::lowest();
        float max_f = std::numeric_limits<float>::max();
        pMin = Vector3f(max_f, max_f, max_f);
        pMax = Vector3f(min_f, min_f, min_f);
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
        //返回包围框在 x、y、z 轴上最长的轴的索引（0 表示 x 轴，1 表示 y 轴，2 表示 z 轴）
        Vector3f d = Diagonal();
        if(d.x() > d.y() && d.x() > d.z())
            return 0;
        else if(d.y() > d.z())
            return 1;
        else
            return 2;
    }

    double SurfaceArea() const
    {
        //计算并返回包围框的表面积
        Vector3f d = Diagonal();
        return 2 * (d.x() * d.y() + d.x() * d.z() + d.y() * d.z());
    }

    Vector3f Centroid()
    {
        return 0.5 * pMin + 0.5 * pMax;
    }

    Bounds3 Intersect(const Bounds3& b)
    {
        Vector3f p1 = Vector3f(fmax(pMin.x(), b.pMin.x()), fmax(pMin.y(), b.pMin.y()), fmax(pMin.z(), b.pMin.z()));
        Vector3f p2 = Vector3f(fmin(pMax.x(), b.pMax.x()), fmin(pMax.y(), b.pMax.y()), fmin(pMax.z(), b.pMax.z()));
        return Bounds3(p1, p2);
    }

    Vector3f Offset(const Vector3f& p) const
    {
        //计算点 p 在包围框内的相对位置，即点 p 在包围框中的归一化坐标
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
        //检查两个包围框是否有重叠部分
        bool x = (b1.pMax.x() >= b2.pMin.x()) && (b1.pMin.x() <= b2.pMax.x());
        bool y = (b1.pMax.y() >= b2.pMin.y()) && (b1.pMin.y() <= b2.pMax.y());
        bool z = (b1.pMax.z() >= b2.pMin.z()) && (b1.pMin.z() <= b2.pMax.z());
        return (x && y && z);
    }

    bool Inside(const Vector3f& p, const Bounds3& b)
    {
        //检查点 p 是否在包围框内
        bool x = p.x() >= b.pMin.x() && p.y() >= pMin.y() && p.z() >= pMin.z();
        bool y = p.x() <= b.pMax.x() && p.y() <= pMax.y() && p.z() <= pMax.z();
        return x && y;
    }

    inline const Vector3f& operator[] (int index) const
    {
        return (index == 0) ? pMin : pMax;
    }

    inline bool IntersectP(const Ray& ray, const Vector3f& invDir,
                           const std::array<int, 3>& dirisNeg) const;
};

inline bool Bounds3::IntersectP(const Ray& ray, const Vector3f& invDir,
                    const std::array<int, 3>& dirisNeg) const
{
    // invDir: ray direction(x,y,z), invDir=(1.0/x,1.0/y,1.0/z), use this because Multiply is faster that Division
    // dirIsNeg: ray direction(x,y,z), dirIsNeg=[int(x>0),int(y>0),int(z>0)], use this to simplify  logic
     //检查射线是否与包围框相交

    double tEnter = -std::numeric_limits<double>::infinity();
    double tExit = std::numeric_limits<double>::infinity();
    for(int i = 0; i < 3; i++)
    {
        double tmin = (pMin[i] - ray.getOrigin()[i]) * invDir[i];
        double tmax = (pMax[i] - ray.getOrigin()[i]) * invDir[i];

        if(!dirisNeg[i])
            std::swap(tmin, tmax);
        // 更新 tenter，取最大值，因为我们需要射线在三个维度上同时进入包围框的最晚时间。
        tEnter = std::max(tEnter, tmin);
        // 更新 texit，取最小值，因为我们需要射线在三个维度上存在离开包围框的最早时间。
        tExit = std::min(tExit, tmax);
    }
    // tenter < texit: 射线进入时间必须早于射线退出时间，这表示射线确实穿过了包围框。
    // texit >= 0: 射线退出时间必须大于等于零，这表示包围框在射线的正向路径上。
    return tEnter <= tExit && tExit >= 0;
}

inline Bounds3 Union(const Bounds3& b1, const Bounds3& b2)
{
    //计算两个包围框的并集。
    Bounds3 res;
    res.pMin = Vector3f::Min(b1.pMin, b2.pMin);
    res.pMax = Vector3f::Max(b1.pMax, b2.pMax);
    return res;
}

inline Bounds3 Union(const Bounds3& b, const Vector3f& p)
{
    //计算一个包围框和一个点的并集。
    Bounds3 res;
    res.pMin = Vector3f::Min(b.pMin, p);
    res.pMax = Vector3f::Max(b.pMax, p);
    return res;
}