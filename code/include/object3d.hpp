#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "ray.hpp"
#include "Material.hpp"
#include "Bounds3.hpp"
#include "global.hpp"
#include "Intersection.hpp"

// Base class for all 3d entities.
class Object3D {
public:
    Material *m;

    Object3D() : m(nullptr) {}

    Object3D(Material* material) : m(material) {}
    virtual ~Object3D() = default;

    // Intersect Ray with this object. If hit, store information in hit structure.
    virtual Intersection getIntersection(Ray _ray) = 0;
    virtual void getSurfaceProperties(const Vector3f &, const Vector3f &, const uint32_t &, const Vector2f &, Vector3f &, Vector2f &) const = 0;
    virtual Vector3f evalDiffuseColor(const Vector2f &) const =0;
    virtual Bounds3 getBounds()=0;
    //与6相比加了新的属性：area，以实现对光源按面积采样
    virtual float getArea() = 0;
    virtual void Sample(Intersection &pos, float &pdf) = 0;
    virtual bool hasEmit() = 0;
    
};

#endif

// #pragma once
// #ifndef RAYTRACING_OBJECT_H
// #define RAYTRACING_OBJECT_H

// #include "Vector3f.h"
// #include "global.hpp"
// #include "Bounds3.hpp"
// #include "ray.hpp"
// #include "Intersection.hpp"

// class Object3D
// {
// public:
//     Object3D() {}
//     virtual ~Object3D() {}
//     //virtual bool intersect(const Ray& ray) = 0;
//     //virtual bool intersect(const Ray& ray, float &, uint32_t &) const = 0;
//     virtual Intersection getIntersection(Ray _ray) = 0;
//     virtual void getSurfaceProperties(const Vector3f &, const Vector3f &, const uint32_t &, const Vector2f &, Vector3f &, Vector2f &) const = 0;
//     //virtual Vector3f evalDiffuseColor(const Vector2f &) const =0;
//     virtual Bounds3 getBounds()=0;
//     //与6相比加了新的属性：area，以实现对光源按面积采样
//     virtual float getArea()=0;
//     virtual void Sample(Intersection &pos, float &pdf)=0;
//     virtual bool hasEmit()=0;
// };



// #endif //RAYTRACING_OBJECT_H
