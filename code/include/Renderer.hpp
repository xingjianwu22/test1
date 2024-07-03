#include "Scene.hpp"
#include "vecmath.h"

#pragma once
struct hit_payload
{
    float tNear;
    uint32_t index;
    Vector2f uv;
    Object3D* hit_obj;
};

class Renderer
{
public:
    void Render(const Scene& scene);

private:
};

