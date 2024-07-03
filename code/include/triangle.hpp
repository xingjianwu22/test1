#pragma once

#include "BVH.hpp"
#include "Intersection.hpp"
#include "Material.hpp"
#include "object3d.hpp"
#include "OBJ_Loader.hpp"
#include "Vector3f.h"
#include <cassert>
#include <array>

bool rayTriangleIntersect(const Vector3f& v0, const Vector3f& v1,
                          const Vector3f& v2, const Vector3f& orig,
                          const Vector3f& dir, float& tnear, float& u, float& v)
{
    Vector3f edge1 = v1 - v0;
    Vector3f edge2 = v2 - v0;
    Vector3f pvec = Vector3f::cross(dir, edge2); //S1
    float det = Vector3f::dot(edge1, pvec);      //S1 * E1
    if (det == 0 || det < 0)
        return false;

    Vector3f tvec = orig - v0;                //S
    u = Vector3f::dot(tvec, pvec);               //S * S1
    if (u < 0 || u > det)
        return false;

    Vector3f qvec = Vector3f::cross(tvec, edge1);//S2
    v = Vector3f::dot(dir, qvec);                //S2 * D
    if (v < 0 || u + v > det)
        return false;

    float invDet = 1 / det;

    tnear = Vector3f::dot(edge2, qvec) * invDet;
    u *= invDet;
    v *= invDet;

    return true;
}

class Triangle : public Object3D
{
public:
    Vector3f v0, v1, v2; // vertices A, B ,C , counter-clockwise order
    Vector3f e1, e2;     // 2 edges v1-v0, v2-v0;
    Vector3f t0, t1, t2; // texture coords
    Vector3f normal;
    Material* m;
    float area;

    // 构造函数，初始化三角形的顶点和材质
    Triangle(Vector3f _v0, Vector3f _v1, Vector3f _v2, Material* _m = nullptr)
        : v0(_v0), v1(_v1), v2(_v2), m(_m)
    {
        e1 = v1 - v0;
        e2 = v2 - v0;
        normal = Vector3f::cross(e1, e2).normalized();
        area = Vector3f::cross(e1, e2).length() * 0.5f;
    }
    // 获取射线与三角形的交点信息
    Intersection getIntersection(Ray ray) override;
    // 获取三角形表面属性
    void getSurfaceProperties(const Vector3f& P, const Vector3f& I,
                              const uint32_t& index, const Vector2f& uv,
                              Vector3f& N, Vector2f& st) const override
    {
        N = normal;
        //        throw std::runtime_error("triangle::getSurfaceProperties not
        //        implemented.");
    }
     // 计算并返回给定纹理坐标处的漫反射颜色
    //Vector3f evalDiffuseColor(const Vector2f&) const override;
    // 获取三角形的包围盒
    Bounds3 getBounds() override;
    void Sample(Intersection &pos, float &pdf){
        float x = std::sqrt(get_random_float()), y = get_random_float();
        pos.coords = v0 * (1.0f - x) + v1 * (x * (1.0f - y)) + v2 * (x * y);
        pos.normal = this->normal;
        pdf = 1.0f / area;
    }
    float getArea(){
        return area;
    }
    bool hasEmit(){
        return m->hasEmission();
    }
    Vector3f evalDiffuseColor(const Vector2f &) const;
};

class MeshTriangle : public Object3D
{
public:
    MeshTriangle(const std::string& filename, Material *mt = new Material())
    {
        objl::Loader loader;
        loader.LoadFile(filename);
        area = 0;
        m = mt;
        std::cout<<"size = "<<loader.LoadedMeshes.size()<<std::endl;
        assert(loader.LoadedMeshes.size() == 1);
        auto mesh = loader.LoadedMeshes[0];

        Vector3f min_vert = Vector3f{std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity()};
        Vector3f max_vert = Vector3f{-std::numeric_limits<float>::infinity(),
                                     -std::numeric_limits<float>::infinity(),
                                     -std::numeric_limits<float>::infinity()};
        for (int i = 0; i < mesh.Vertices.size(); i += 3) {
            std::array<Vector3f, 3> face_vertices;

            for (int j = 0; j < 3; j++) {
                auto vert = Vector3f(mesh.Vertices[i + j].Position.X,
                                     mesh.Vertices[i + j].Position.Y,
                                     mesh.Vertices[i + j].Position.Z);
                face_vertices[j] = vert;

                min_vert = Vector3f(std::min(min_vert.x(), vert.x()),
                                    std::min(min_vert.y(), vert.y()),
                                    std::min(min_vert.z(), vert.z()));
                max_vert = Vector3f(std::max(max_vert.x(), vert.x()),
                                    std::max(max_vert.y(), vert.y()),
                                    std::max(max_vert.z(), vert.z()));
            }

            triangles.emplace_back(face_vertices[0], face_vertices[1],
                                   face_vertices[2], mt);
        }

        bounding_box = Bounds3(min_vert, max_vert);

        std::vector<Object3D*> ptrs;
        for (auto& tri : triangles){
            ptrs.push_back(&tri);
            area += tri.area;
        }
        bvh = new BVHAccel(ptrs);
    }


    Bounds3 getBounds() { return bounding_box; }

    void getSurfaceProperties(const Vector3f& P, const Vector3f& I,
                              const uint32_t& index, const Vector2f& uv,
                              Vector3f& N, Vector2f& st) const
    {
        const Vector3f& v0 = vertices[vertexIndex[index * 3]];
        const Vector3f& v1 = vertices[vertexIndex[index * 3 + 1]];
        const Vector3f& v2 = vertices[vertexIndex[index * 3 + 2]];
        Vector3f e0 = (v1 - v0).normalized();
        Vector3f e1 = (v2 - v1).normalized();
        N = Vector3f::cross(e0, e1).normalized();
        const Vector2f& st0 = stCoordinates[vertexIndex[index * 3]];
        const Vector2f& st1 = stCoordinates[vertexIndex[index * 3 + 1]];
        const Vector2f& st2 = stCoordinates[vertexIndex[index * 3 + 2]];
        st = st0 * (1 - uv.x() - uv.y()) + st1 * uv.x() + st2 * uv.y();
    }

    Vector3f evalDiffuseColor(const Vector2f& st) const
    {
        float scale = 5;
        float pattern =
            (fmodf(st.x() * scale, 1) > 0.5) ^ (fmodf(st.y() * scale, 1) > 0.5);
        Vector3f v1 = Vector3f(0.815, 0.235, 0.031);
        Vector3f v2 = Vector3f(0.937, 0.937, 0.231);
        return v1 * (1 - pattern) + v2 * pattern;
    }

    Intersection getIntersection(Ray ray)
    {
        Intersection intersec;

        if (bvh) {
            intersec = bvh->Intersect(ray);
        }

        return intersec;
    }

    void Sample(Intersection &pos, float &pdf){
        bvh->Sample(pos, pdf);
        pos.emit = m->getEmission();
    }
    float getArea(){
        return area;
    }
    bool hasEmit(){
        return m->hasEmission();
    }


    Bounds3 bounding_box;
    std::unique_ptr<Vector3f[]> vertices;
    uint32_t numTriangles;
    std::unique_ptr<uint32_t[]> vertexIndex;
    std::unique_ptr<Vector2f[]> stCoordinates;

    std::vector<Triangle> triangles;

    BVHAccel* bvh;
    float area;

    Material* m;
};


inline Bounds3 Triangle::getBounds() { return Union(Bounds3(v0, v1), v2); }

inline Intersection Triangle::getIntersection(Ray ray)
{
    Intersection inter;

    if (Vector3f::dot(ray.getDirection(), normal) > 0)
        return inter;
    double u, v, t_tmp = 0;
    Vector3f pvec = Vector3f::cross(ray.getDirection(), e2);  //S1
    double det = Vector3f::dot(e1, pvec);                //S1 * E1
    if (fabs(det) < EPSILON)                          //det趋近于0，光线与面平行
        return inter;

    double det_inv = 1. / det;                        //1 / DET
    Vector3f tvec = ray.getOrigin() - v0;                  //S
    u = Vector3f::dot(tvec, pvec) * det_inv;             //S * S1 / DET
    if (u < 0 || u > 1)
        return inter;
    Vector3f qvec = Vector3f::cross(tvec, e1);           //S2 = S X E1
    v = Vector3f::dot(ray.getDirection(), qvec) * det_inv;    //S2 * D
    if (v < 0 || u + v > 1)
        return inter;
    t_tmp = Vector3f::dot(e2, qvec) * det_inv;           //E2 * S2 / DET
    if(t_tmp < 0)
        return inter;
    inter.isIntersected = true;
    inter.coords = ray(t_tmp);
    inter.normal = normal;
    inter.distance = t_tmp;
    inter.obj = this;
    inter.m = m;
    return inter;


    return inter;
}

inline Vector3f Triangle::evalDiffuseColor(const Vector2f&) const
{
    return Vector3f(0.5, 0.5, 0.5);
}

