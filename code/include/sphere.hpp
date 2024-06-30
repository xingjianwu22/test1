#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include "Bounds3.hpp"
#include "Material.hpp"
#include <vecmath.h>
#include <cmath>

// TODO: Implement functions and add more fields as necessary

class Sphere : public Object3D {
public:
    Vector3f center;
    float radius, radius2;
    float area;
    Sphere() : radius(1), radius2(1), center(0, 0, 0) {
        // unit ball at the center
    }

    Sphere(const Vector3f &center, float radius, Material *material)
    {
        this->center = center;
        this->radius = radius;
        this->radius2 = radius * radius;
        this->m = material;
        this->area = 4 * M_PI * radius2;
    }

    ~Sphere() override = default;

    Intersection getIntersection(Ray ray){
        Intersection result;
        result.isIntersected = false;
        Vector3f L = ray.getOrigin() - center;
        float a = Vector3f::dot(ray.getDirection(), ray.getDirection());
        float b = 2 * Vector3f::dot(ray.getDirection(), L);
        float c = Vector3f::dot(L, L) - radius2;
        float t0, t1;
        if (!solveQuadratic(a, b, c, t0, t1)) return result;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) 
            return result;
        
        result.isIntersected=true;
        result.coords = ray.pointAtParameter(t0);
        result.normal = Vector3f(result.coords - center).normalized();
        result.m = this->m;
        result.obj = this;
        result.distance = t0;
        return result;
    }
    void getSurfaceProperties(const Vector3f &P, const Vector3f &I, const uint32_t &index, const Vector2f &uv, Vector3f &N, Vector2f &st) const
    { 
        N = (P - center).normalized(); 
    }

    Vector3f evalDiffuseColor(const Vector2f &st)const {
        return m->getDiffuseColor();
    }

    Bounds3 getBounds(){
        return Bounds3(Vector3f(center.x() - radius, center.y() - radius, center.z() - radius),
                       Vector3f(center.x() + radius, center.y() + radius, center.z() + radius));
    }

    // 将光源进行按面采样,随机从光源发射一条ray打到场景中的sphere上得到某个交点
    void Sample(Intersection &pos, float &pdf){
        float theta = 2.0 * M_PI * get_random_float(), phi = M_PI * get_random_float();
        Vector3f dir(std::cos(phi), std::sin(phi)*std::cos(theta), std::sin(phi)*std::sin(theta));
        pos.coords = center + radius * dir;
        pos.normal = dir;
        pos.emit = m->getEmission();
        pdf = 1.0f / area;
    }
    float getArea(){
        return area;
    }
    bool hasEmit(){
        return m->hasEmission();
    }
};

#endif
