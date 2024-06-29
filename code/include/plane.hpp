#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement Plane representing an infinite plane
// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() : normal(Vector3f::UP), d(0) {}

    Plane(const Vector3f &normal, float d, Material *m)
    {
        this -> normal = normal;
        this -> d = d;
        this -> m = m;
    }

    ~Plane() override = default;

    Intersection getIntersection(Ray ray)
    {
        Intersection result;
        result.isIntersected = false;
        Vector3f o(ray.getOrigin()), dir(ray.getDirection());
        float cos = Vector3f::dot(normal, dir);
        //平行
        if (fabs(cos) < 1e-6)
            return result;
        float t = (d - Vector3f::dot(o, normal)) / Vector3f::dot(normal, dir);
        if(t < 0)
            return result;
        else
        {
            result.isIntersected = true;
            result.coords = ray.pointAtParameter(t);
            result.m = this->m;
            result.obj = this;
            result.distance = t;
        }
        return result;
    }

protected:
    Vector3f normal; // (a, b, c)
    float d;         // d
    //平面方程：ax + by + cz = d;
};

#endif //PLANE_H
		

