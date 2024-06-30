
#pragma once
#include "Vector3f.h"
#include "Material.hpp"
class Object3D;
class Sphere;

struct Intersection
{
    Intersection(){
        isIntersected = false;
        coords=Vector3f();
        normal=Vector3f();
        distance= std::numeric_limits<double>::max();
        obj = nullptr;
        m = nullptr;
    }
    bool isIntersected;
    Vector3f coords;
    Vector3f tcoords;
    Vector3f normal;
    Vector3f emit;
    double distance;
    Object3D* obj;
    Material* m;
};
