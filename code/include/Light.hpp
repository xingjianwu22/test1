#pragma once

#include "Vector3f.h"
#include "global.hpp"

class Light
{
public:
    Light(const Vector3f &p, const Vector3f &i) : position(p), intensity(i) {}
    virtual ~Light() = default;
    Vector3f position;
    Vector3f intensity;
};


// AreadLight 用于表示一个有一定面积的区域光源
class AreaLight : public Light
{
public:
    float length;
    Vector3f normal;
    Vector3f u;
    Vector3f v;

    AreaLight(const Vector3f &p, const Vector3f &i) : Light(p, i)
    {
        normal = Vector3f(0, -1, 0);
        u = Vector3f(1, 0, 0);   //u和v表示正方形光源两条边的方向向量
        v = Vector3f(0, 0, 1);
        length = 100;            //边的长度
    }

    Vector3f SamplePoint() const
    {
        auto random_u = get_random_float();
        auto random_v = get_random_float();
        return position + random_u * u + random_v * v;
    }

};




