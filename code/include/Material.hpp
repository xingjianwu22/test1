#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "global.hpp"
#include <iostream>

// TODO: Implement Shade function that computes Phong introduced in class.

enum MaterialType { DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION, REFRACTION, DIFFUSE, GLOSSY, MICROFACET };

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
    float roughness;
    Vector3f albedo;
    Vector3f F0;
    float metallic;

    // explicit Material(MaterialType t = DIFFUSE_AND_GLOSSY, const Vector3f &d_color = Vector3f::ZERO, const Vector3f &s_color = Vector3f::ZERO, float s = 0) :
    //         m_type(t), diffuseColor(d_color), specularColor(s_color), shininess(s) {

    // }
    Material(MaterialType t = DIFFUSE, Vector3f e = Vector3f::ZERO) : m_type(t), m_emission(e) {}
    virtual ~Material() = default;

    Material(MaterialType t, const Vector3f &e, const Vector3f &albedo, float roughness, float metallic) {
        m_type = t;
        m_emission = e;
        this->albedo = albedo;
        this->roughness = roughness;
        this->metallic = metallic;
        Vector3f base(0.04);
        F0 = Vector3f::lerp(base, albedo, metallic);
    }
    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }
    MaterialType getType() {return m_type;}
    Vector3f getEmission() {return m_emission;}
    bool hasEmission()
    {
        if(m_emission.length() < 0.00001)
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
    Vector3f reflect(const Vector3f &I, const Vector3f &N) const
    {
        return I - 2 * Vector3f::dot(I, N) * N;
    }

    Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior) const
    {
        float cosi = clamp(-1, 1, Vector3f::dot(I, N));
        float etai = 1, etat = ior;
        Vector3f n = N;
        if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
    }

    void fresnel(const Vector3f &I, const Vector3f &N, const float &ior, float &kr) const
    {
        float cosi = clamp(-1, 1, Vector3f::dot(I, N));
        float etai = 1, etat = ior;
        if (cosi > 0) {  std::swap(etai, etat); }
        // Compute sini using Snell's law
        float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
        // Total internal reflection
        if (sint >= 1) {
            kr = 1;
        }
        else {
            float cost = sqrtf(std::max(0.f, 1 - sint * sint));
            cosi = fabsf(cosi);
            float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
            float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
            kr = (Rs * Rs + Rp * Rp) / 2;
        }
        // As a consequence of the conservation of energy, transmittance is given by:
        // kt = 1 - kr;
    }

    Vector3f importanceSampling(const Vector3f &wo, const Vector3f &normal) {
        double r0 = get_random_float();
        double r1 = get_random_float();
        double a2 = roughness * roughness;
        double theta = acos(sqrt((1 - r0) / ((a2 - 1) * r0 + 1)));
        double phi = 2 * M_PI * r1;
        
        // double x = sin(theta) * cos(phi);
        // double y = cos(theta);
        // double z = sin(theta) * sin(phi);
        double x = sin(theta) * cos(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(theta);
        Vector3f wm = Vector3f(x, y, z);
        Vector3f wm_w = toWorld(wm, normal);
        return reflect(wo, wm_w);
    }
    
    double importancePDF(const Vector3f &wo, const Vector3f &wi, const Vector3f &normal) {
        Vector3f h = (wo + wi).normalized();
        double cosTheta = Vector3f::dot(normal, h);
        double D = distributionGGX(normal, h, roughness);
        return (D * cosTheta) / (4.0f * Vector3f::dot(wo, h));
    }
    // sample a ray by Material properties
    inline Vector3f sample(const Vector3f &wi, const Vector3f &N);
    // given a ray, calculate the PdF of this ray
    inline float pdf(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    // given a ray, calculate the contribution of this ray
    inline Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N);
    
private:
    inline double distributionGGX(const Vector3f &normal, const Vector3f &h, double rough_ness) {
        double a2 = rough_ness * rough_ness;
        double nDotH = std::max(Vector3f::dot(normal, h), 0.0f);
        double nDotH2 = nDotH * nDotH;

        double denom = nDotH2 * (a2 - 1.0f) + 1.0f;
        denom = M_PI * denom * denom;
        return a2 / denom;
    }

    inline double geometrySchlickGGX(double nDotV, double k) {
        double denom = nDotV * (1.0f - k) + k;
        return nDotV / denom;
    }

    inline double geometrySmith(const Vector3f &normal, const Vector3f &v, const Vector3f &l, double k) {
        double nDotV = std::max(Vector3f::dot(normal, v), 0.0f);
        double nDotL = std::max(Vector3f::dot(normal, l), 0.0f);
        double ggx1 = geometrySchlickGGX(nDotV, k);
        double ggx2 = geometrySchlickGGX(nDotL, k);

        return ggx1 * ggx2;
    }

    inline Vector3f fresnelSchlick(double cosTheta, const Vector3f &F0) {
        return F0 + (Vector3f(1.0f) - F0) * pow(1.0 - cosTheta, 5.0);
    }
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
        case MICROFACET:
        {
            return importanceSampling(wi, N);
        }
        default:
            return Vector3f(0.0f);
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
        case MICROFACET:
        {
            return std::max((double)0.0001f, importancePDF(wo, wi, N));
        }
        default:
            return 0.0f;
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
        case REFLECTION: {
            Vector3f reflectDir = reflect(-wi, N);
            if (Vector3f::dot(reflectDir, wo) > 0.0f) {
                return Ks / Vector3f::dot(N, wo);
            }
            return Vector3f(0.0f);
        }
        case REFRACTION: {
            float F;
            fresnel(wi, N, ior, F);
            Vector3f refractDir = refract(-wi, N, ior);
            if (Vector3f::dot(refractDir, wo) > 0.0f) {
                return (1.0f - F) * Ks / Vector3f::dot(N, wo);
            }
            return Vector3f(0.0f);
        }
        case GLOSSY:
        case MICROFACET:
        {
            double cos1 = std::max(Vector3f::dot(N, wo), 0.0f);
            double cos2 = std::max(Vector3f::dot(N, wi), 0.0f);
            if (cos1 > 0.0f && cos2 > 0.0f) {
                Vector3f h = (wi + wo).normalized();
                double k = pow((roughness + 1.0f), 2) / 8.0f;
                double distribute = distributionGGX(N, h, roughness);
                double geometry = geometrySmith(N, wo, wi, k);
                
                Vector3f fresnel = fresnelSchlick(cos1, F0);
                Vector3f Ks = fresnel;
                Vector3f Kd = (Vector3f(1.0f) - Ks) * (1.0f - metallic);
                return Kd * albedo / M_PI + Ks * distribute * geometry / std::max((double)0.0001f, (4.0f * cos1 * cos2));
            }
            else
                return Vector3f(0.0f);
            break;
        }
        default:
            return Vector3f(0.0f);
    }
}



#endif // MATERIAL_H



