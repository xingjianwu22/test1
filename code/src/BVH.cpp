#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object3D*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    //构造函数，用于初始化 BVHAccel 对象并构建 BVH 树。
    time_t start, stop;
    time(&start); //获取当前时间并将其存储在 start 变量中。
    //具体来说，time 函数返回从1970年1月1日00:00:00 UTC到当前时间的秒数，即所谓的 "Unix 时间戳"
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object3D*> objects)
{
    //递归地构建 BVH 树，将对象组织成层次结构。
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        node->area = objects[0]->getArea();
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector<Object3D*>{objects[0]});
        node->right = recursiveBuild(std::vector<Object3D*>{objects[1]});
        node->area = node->left->area + node->right->area;
        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x() <
                       f2->getBounds().Centroid().x();
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y() <
                       f2->getBounds().Centroid().y();
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z() <
                       f2->getBounds().Centroid().z();
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object3D*>(beginning, middling);
        auto rightshapes = std::vector<Object3D*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
        node->area = node->left->area + node->right->area;
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    //检测射线 ray 是否与 BVH 中的任何对象相交，并返回最近的相交点信息。
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    //递归地检测射线 ray 是否与当前 BVH 树节点 node 及其子节点中的任何对象相交，并返回最近的相交点信息。
    
    // TODO Traverse the BVH to find intersection
    //递归调用Bounds3::intersectP
    Intersection res;
    // 首先判断这个ray与当前box是否有交点：
    // 1.如果没有交点 -> 那就不用继续了，因为再进行细分也没有意义，直接返回当前的intersection
    std::array<int, 3> dirIsNeg = { int(ray.getDirection().x() > 0),int(ray.getDirection().y() > 0),int(ray.getDirection().z() > 0) };
    if (!node->bounds.IntersectP(ray, ray.getInvDirection(), dirIsNeg)) {
        return res;
    }
    // 2.如果有交点 -> 2.1若该点无子节点，则返回该交点，该节点是没有子节点，可以进行BVN判断是否相交
    if (node->left == nullptr && node->right == nullptr) {
        res = node->object->getIntersection(ray);
        return res;
    }
    //2.1该点有子节点，则左右子节点分别判断，继续递归
    Intersection resleft, resright;
    resleft = getIntersection(node->left, ray);
    resright = getIntersection(node->right, ray);
    return resleft.distance < resright.distance ? resleft : resright;
}

void BVHAccel::getSample(BVHBuildNode* node, float p, Intersection &pos, float &pdf){
    if(node->left == nullptr || node->right == nullptr){
        node->object->Sample(pos, pdf);
        pdf *= node->area;
        return;
    }
    if(p < node->left->area)
        getSample(node->left, p, pos, pdf);
    else
        getSample(node->right, p - node->left->area, pos, pdf);
}

void BVHAccel::Sample(Intersection &pos, float &pdf){
    float p = std::sqrt(get_random_float()) * root->area;
    getSample(root, p, pos, pdf);
    pdf /= root->area;
}