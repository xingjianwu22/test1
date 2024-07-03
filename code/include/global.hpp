#pragma once
#include <iostream>
#include <cmath>
#include <random>

#undef M_PI  //#undef M_PI 确保之前可能在其他头文件中定义的 M_PI 宏被取消定义。
#define M_PI 3.141592653589793f //重新定义

extern const float  EPSILON;
const float kInfinity = std::numeric_limits<float>::max();

inline float clamp(const float &lo, const float &hi, const float &v)
{ return std::max(lo, std::min(hi, v)); }

inline  bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0) return false;
    else if (discr == 0) x0 = x1 = - 0.5 * b / a;
    else {
        float q = (b > 0) ?
                  -0.5 * (b + sqrt(discr)) :
                  -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1) std::swap(x0, x1);
    return true;
}

inline float get_random_float()
{
    // std::random_device 是一个随机数生成器，通常用于获取硬件生成的随机数种子。
    static std::random_device dev;
    // std::mt19937 是一种梅森旋转算法的伪随机数生成器，使用 dev() 生成的随机种子进行初始化。
    static std::mt19937 rng(dev());
    // std::uniform_real_distribution<float> 是一个均匀分布的随机数分布器
    static std::uniform_real_distribution<float> dist(0.f, 1.f); // distribution in range [0, 1)
    return dist(rng);
}

inline void UpdateProgress(float progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
};

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
