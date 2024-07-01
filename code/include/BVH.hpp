#pragma once

#include <atomic>
#include <vector>
#include <memory>
#include <ctime>
#include "object3d.hpp"
#include "ray.hpp"
#include "Bounds3.hpp"
#include "Intersection.hpp"
#include "Vector3f.h"

struct BVHBuildNode;
// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;

// BVHAccel Declarations
inline int leafNodes, totalLeafNodes, totalPrimitives, interiorNodes;
class BVHAccel {

public:
    // BVHAccel Public Types
    enum class SplitMethod { NAIVE, SAH }; //定义了分割方法，NAIVE 表示简单的分割方法，SAH 表示表面面积启发式（Surface Area Heuristic）。

    // BVHAccel Public Methods
    //maxPrimsInNode: 每个叶节点中最多包含的对象数。
    BVHAccel(std::vector<Object3D*> p, int maxPrimsInNode = 1, SplitMethod splitMethod = SplitMethod::NAIVE);
    //返回整个 BVH 的包围盒。
    Bounds3 WorldBound() const;
    ~BVHAccel();

    //检测射线 ray 是否与 BVH 中的任何对象相交，并返回最近的交点信息
    Intersection Intersect(const Ray &ray) const;
    //递归检测射线 ray 是否与 BVH 子树中的任何对象相交，并返回最近的交点信息。
    Intersection getIntersection(BVHBuildNode* node, const Ray& ray)const;
    //BVH 树的根节点
    BVHBuildNode* root;

    // BVHAccel Private Methods
    BVHBuildNode* recursiveBuild(std::vector<Object3D*>objects);

    // BVHAccel Private Data
    //每个叶节点中最多包含的对象数。
    const int maxPrimsInNode;
    //分割策略
    const SplitMethod splitMethod;
    //场景中的对象列表。
    std::vector<Object3D*> primitives;

    void getSample(BVHBuildNode* node, float p, Intersection &pos, float &pdf);
    void Sample(Intersection &pos, float &pdf);
};

struct BVHBuildNode {
    Bounds3 bounds;        // 节点的包围盒
    BVHBuildNode *left;    // 左子节点指针
    BVHBuildNode *right;   // 右子节点指针
    Object3D* object;      // 节点包含的对象（仅在叶节点中使用）
    float area;            // 面积
 
public:
    int splitAxis=0;       //分割轴
    int firstPrimOffset=0; //第一个对象的偏移
    int nPrimitives=0;     //节点中包含的对象数
    // BVHBuildNode Public Methods
    BVHBuildNode(){
        bounds = Bounds3();
        left = nullptr;right = nullptr;
        object = nullptr;
    }
};




