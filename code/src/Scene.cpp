#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    //1.判断是否有交点：光线与场景中物体相交？
    Intersection inter = Scene::intersect(ray);
    //如果没交点，返回背景光
    if (inter.isIntersected)
        return shade(inter, ray);

    return this->backgroundColor;
}

Vector3f Scene::shade(const Intersection &isec, const Ray &ray) const
{
    // TO DO Implement Path Tracing Algorithm here

    //创建变量以储存直接和间接光照计算值
    Vector3f dir = { 0.0,0.0,0.0 };
    Vector3f indir = { 0.0,0.0,0.0 };
    Material *m = isec.m;
    //2.ray打到光源了：说明渲染方程只用算前面的自发光项，因此直接返回材质的自发光项
    if (m->hasEmission())
        return m->getEmission();

    //3.ray打到物体
    //对场景中的光源进行采样，得到采样点light_pos和pdf_light
    Intersection light_pos;
    float pdf_light = 0.0f;
    sampleLight(light_pos, pdf_light);
   
    //3.1计算直接光照
 
    //物体的一些参数
    Vector3f p = isec.coords;
    Vector3f N = isec.normal.normalized();
    Vector3f wo = ray.getDirection();//物体指向场景
    //光源的一些参数
    Vector3f xx = light_pos.coords;
    Vector3f NN = light_pos.normal.normalized();
    Vector3f ws = (xx - p).normalized();//物体指向光源
    float dis = (p - xx).length();//二者距离
    float dis2 = Vector3f::dot((p - xx), (p - xx));
    
    //判断光源与物体间是否有遮挡：
    //发出一条射线，方向为ws 物体p -> 光源xx
    Vector3f p_deviation = (Vector3f::dot(ray.getDirection(), N) < 0) ?
        p + N * EPSILON :
        p - N * EPSILON ;
    Ray light_to_obj(p_deviation, ws);//Ray(orig,dir)
    Intersection light_to_scene = Scene::intersect(light_to_obj);
    bool hasDirectLight = false;
    //假如dis>light_to_scene.distance就说明有遮挡，那么反着给条件即可：
    if (light_to_scene.isIntersected&& (light_to_scene.distance-dis>-sqrt(EPSILON))) {//没有遮挡
        // 计算直接光照
        Vector3f L_i = light_pos.emit;//光强
        Vector3f f_r = isec.m->eval(ws, -wo, N);//BRDF==材质
        float cos_theta = Vector3f::dot(ws, N);//物体夹角
        float cos_theta_l = Vector3f::dot(-ws, NN);//光源夹角
        dir = L_i * f_r * cos_theta * cos_theta_l / dis2 / pdf_light;
        // due to sometimes we won't get the correct ray to reflect light in the specular surface.
        // thus, we let the ray to choose it's own directon.
        hasDirectLight = dir.length() > 0.1; 
    }
    //3.2间接光照
    
    //俄罗斯轮盘赌
    //Scene.hpp中已经定义了P_RR:RussianRoulette=0.8
    float ksi = get_random_float();//随机取[0,1]
    if (ksi < RussianRoulette) {
        //计算间接光照
        
        //随机生成一个wi方向
        Vector3f wi = isec.m->sample(wo, N).normalized();//这里的wi其实没参与计算，返回的是一个随机的方向
        Ray r(p, wi);
        Intersection obj_to_scene = Scene::intersect(r);
        //击中了物体 && 物体不是光源
        if (obj_to_scene.isIntersected && (!obj_to_scene.m->hasEmission() || !hasDirectLight)) {
            Vector3f f_r = m->eval(wi, -wo, N);
            float cos_theta = std::max(0.0f, Vector3f::dot(wi, N));
            float pdf_hemi = m->pdf(wi, -wo, N);
            if(pdf_hemi > EPSILON) //防止pdf_hemi取接近0，产生白色噪点
            indir = shade(obj_to_scene, r) * f_r * cos_theta / pdf_hemi / RussianRoulette;
        }
    }
    return Vector3f::Max(Vector3f::Min(dir + indir, Vector3f(1.0f)), Vector3f(0.0f));
}


// Implementation of the Whitted-syle light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected object is either reflective or reflective and refractive,
// then we compute the reflection/refracton direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refractin depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is duffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.
// Vector3f Scene::castRay(const Ray &ray, int depth) const
// {
//     if (depth > this->maxDepth) {
//         return Vector3f(0.0,0.0,0.0);
//     }
//     Intersection intersection = Scene::intersect(ray);
//     Material *m = intersection.m;
//     Object3D *hitObject = intersection.obj;
//     Vector3f hitColor = this->backgroundColor;
// //    float tnear = kInfinity;
//     Vector2f uv;
//     uint32_t index = 0;
//     if(intersection.isIntersected) {

//         Vector3f hitPoint = intersection.coords;
//         Vector3f N = intersection.normal; // normal
//         Vector2f st; // st coordinates
//         hitObject->getSurfaceProperties(hitPoint, ray.getDirection(), index, uv, N, st);
// //        Vector3f tmp = hitPoint;
//         switch (m->getType()) {
//             case REFLECTION_AND_REFRACTION:
//             {
//                 Vector3f reflectionDirection = (reflect(ray.getDirection(), N)).normalized();
//                 Vector3f refractionDirection = (refract(ray.getDirection(), N, m->ior)).normalized();
//                 Vector3f reflectionRayOrig = (Vector3f::dot(reflectionDirection, N) < 0) ?
//                                              hitPoint - N * EPSILON :
//                                              hitPoint + N * EPSILON;
//                 Vector3f refractionRayOrig = (Vector3f::dot(refractionDirection, N) < 0) ?
//                                              hitPoint - N * EPSILON :
//                                              hitPoint + N * EPSILON;
//                 Vector3f reflectionColor = castRay(Ray(reflectionRayOrig, reflectionDirection), depth + 1);
//                 Vector3f refractionColor = castRay(Ray(refractionRayOrig, refractionDirection), depth + 1);
//                 float kr;
//                 fresnel(ray.getDirection(), N, m->ior, kr);
//                 hitColor = reflectionColor * kr + refractionColor * (1 - kr);
//                 break;
//             }
//             case REFLECTION:
//             {
//                 float kr;
//                 fresnel(ray.getDirection(), N, m->ior, kr);
//                 Vector3f reflectionDirection = reflect(ray.getDirection(), N);
//                 Vector3f reflectionRayOrig = (Vector3f::dot(reflectionDirection, N) < 0) ?
//                                              hitPoint + N * EPSILON :
//                                              hitPoint - N * EPSILON;
//                 hitColor = castRay(Ray(reflectionRayOrig, reflectionDirection),depth + 1) * kr;
//                 break;
//             }
//             default:
//             {
//                 // [comment]
//                 // We use the Phong illumation model int the default case. The phong model
//                 // is composed of a diffuse and a specular reflection component.
//                 // [/comment]
//                 Vector3f lightAmt = 0, specularColor = 0;
//                 Vector3f shadowPointOrig = (Vector3f::dot(ray.getDirection(), N) < 0) ?
//                                            hitPoint + N * EPSILON :
//                                            hitPoint - N * EPSILON;
//                 // [comment]
//                 // Loop over all lights in the scene and sum their contribution up
//                 // We also apply the lambert cosine law
//                 // [/comment]
//                 for (uint32_t i = 0; i < get_lights().size(); ++i)
//                 {
//                     auto area_ptr = dynamic_cast<AreaLight*>(this->get_lights()[i].get());
//                     if (area_ptr)
//                     {
//                         // Do nothing for this assignment
//                     }
//                     else
//                     {
//                         Vector3f lightDir = get_lights()[i]->position - hitPoint;
//                         // square of the distance between hitPoint and the light
//                         float lightDistance2 = Vector3f::dot(lightDir, lightDir);
//                         lightDir = lightDir.normalized();
//                         float LdotN = std::max(0.f, Vector3f::dot(lightDir, N));
//                         Object3D *shadowHitObject = nullptr;
//                         float tNearShadow = kInfinity;
//                         // is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
//                         bool inShadow = bvh->Intersect(Ray(shadowPointOrig, lightDir)).isIntersected;
//                         lightAmt += (1 - inShadow) * get_lights()[i]->intensity * LdotN;
//                         Vector3f reflectionDirection = reflect(-lightDir, N);
//                         specularColor += powf(std::max(0.f, -Vector3f::dot(reflectionDirection, ray.getDirection())),
//                                               m->shininess) * get_lights()[i]->intensity;
//                     }
//                 }
//                 hitColor = lightAmt * (hitObject->evalDiffuseColor(st) * m->Kd + specularColor * m->Ks);
//                 break;
//             }
//         }
//     }

//     return hitColor;
// }
