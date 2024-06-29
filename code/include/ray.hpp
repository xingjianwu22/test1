#ifndef RAY_H
#define RAY_H

#include <cassert>
#include <iostream>
#include <limits>
#include <Vector3f.h>


// Ray class mostly copied from Peter Shirley and Keith Morley
class Ray {
public:

    Ray() = delete;
    Ray(const Vector3f &orig, const Vector3f &dir, const double _t = 0.0) {
        origin = orig;
        direction = dir;
        t = _t;
        direction_inv = Vector3f(1./direction.x(), 1./direction.y(), 1./direction.z());
        t_min = 0.0;
        t_max = std::numeric_limits<double>::max();
    }

    Ray(const Ray &r) {
        origin = r.origin;
        direction = r.direction;
    }

    const Vector3f &getOrigin() const {
        return origin;
    }

    const Vector3f &getDirection() const {
        return direction;
    }

    const Vector3f &getInvDirection() const {
        return direction_inv;
    }

    Vector3f pointAtParameter(float t) const {
        return origin + direction * t;
    }

    Vector3f operator() (double t) const {
        return origin + t * direction;
    }

    
private:
    Vector3f origin;
    Vector3f direction, direction_inv;
    double t;
    double t_min, t_max;

};

inline std::ostream &operator<<(std::ostream &os, const Ray &r) {
    os << "Ray <" << r.getOrigin() << ", " << r.getDirection() << ">";
    return os;
}

#endif // RAY_H
