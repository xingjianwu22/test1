#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "global.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.

enum MaterialType { DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION, DIFFUSE };

class Material {
public:
    MaterialType m_type;
    Vector3f m_emission;   // 材料自发光的颜色
    Vector3f diffuseColor; // 漫反射颜色，决定材料在漫反射光下的颜色
    Vector3f specularColor;// 镜面反射颜色
    float shininess;       // 高光系数(影响镜面反射的大小和强度)
    Vector3f Kd;           // 漫反射系数，控制材料的漫反射光的强度
    Vector3f Ks;           // 镜面反射系数，控制材料的镜面反射光的强度
    float ior;             // 折射率

    explicit Material(MaterialType t = DIFFUSE_AND_GLOSSY, const Vector3f &d_color = Vector3f::ZERO, const Vector3f &s_color = Vector3f::ZERO, float s = 0) :
            m_type(t), diffuseColor(d_color), specularColor(s_color), shininess(s) {

    }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }
    MaterialType getType() {return m_type;}
    Vector3f getEmission() {return m_emission;}
    bool hasEmission()
    {
        if(m_emission.length() > 0.00001)
            return false;
        else
            return true;
    }
    Vector3f toWorld(const Vector3f &a, const Vector3f &N){
        Vector3f B, C;
        if (std::fabs(N.x()) > std::fabs(N.y())){
            float invLen = 1.0f / std::sqrt(N.x() * N.x() + N.z() * N.z());
            C = Vector3f(N.z() * invLen, 0.0f, -N.x() *invLen);
        }
        else {
            float invLen = 1.0f / std::sqrt(N.y() * N.y() + N.z() * N.z());
            C = Vector3f(0.0f, N.z() * invLen, -N.y() *invLen);
        }
        B = Vector3f::cross(C, N);
        return a.x() * B + a.y() * C + a.z() * N;
    }
    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
};

Vector3f Material::sample(const Vector3f &wi, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample on the hemisphere
            float x_1 = get_random_float(), x_2 = get_random_float();
            float z = std::fabs(1.0f - 2.0f * x_1);
            float r = std::sqrt(1.0f - z * z), phi = 2 * M_PI * x_2;
            Vector3f localRay(r*std::cos(phi), r*std::sin(phi), z);
            return toWorld(localRay, N);            
            break;
        }
    }
}

float Material::pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // uniform sample probability 1 / (2 * PI)
            if (Vector3f::dot(wo, N) > 0.0f)
                return 0.5f / M_PI;
            else
                return 0.0f;
            break;
        }
    }
}

Vector3f Material::eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
    switch(m_type){
        case DIFFUSE:
        {
            // calculate the contribution of diffuse   model
            float cosalpha = Vector3f::dot(N, wo);
            if (cosalpha > 0.0f) {
                Vector3f diffuse = Kd / M_PI;
                return diffuse;
            }
            else
                return Vector3f(0.0f);
            break;
        }
    }
}

#endif // MATERIAL_H
