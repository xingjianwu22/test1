#pragma once

#include "Vector3f.h"
#include "Material.hpp"

struct Object3D;

struct Intersection
{
public:
    bool isIntersected;
    Vector3f coords;
    Vector3f tcoords;
    Vector3f emit;
    Vector3f normal;
    double distance;
    Object3D* obj;
    Material* m;

    Intersection()
    {
        isIntersected = false;
        coords = Vector3f();
        normal = Vector3f();
        distance = std::numeric_limits<double>::max();
        obj = nullptr;
        m = nullptr;
    }


};