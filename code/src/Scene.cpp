//
// Created by Göksu Güvendiren on 2019-05-14.
//

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

// // Implementation of Path Tracing
// Vector3f Scene::castRay(const Ray& ray, int depth) const
// {
//     // TO DO Implement Path Tracing Algorithm here
//     Vector3f hitColor = this->backgroundColor;
//     Intersection shade_point_inter = Scene::intersect(ray);
//     if (shade_point_inter.isIntersected)
//     {

//         Vector3f p = shade_point_inter.coords;
//         Vector3f wo = ray.getDirection();
//         Vector3f N = shade_point_inter.normal;
//         Vector3f L_dir(0), L_indir(0);

//        //sampleLight(inter,pdf_light)
//         Intersection light_point_inter;
//         float pdf_light;
//         sampleLight(light_point_inter, pdf_light);
//         //Get x,ws,NN,emit from inter
//         Vector3f x = light_point_inter.coords;
//         Vector3f ws = (x-p).normalized();
//         Vector3f NN = light_point_inter.normal;
//         Vector3f emit = light_point_inter.emit;
//         float distance_pTox = (x - p).length();
//         //Shoot a ray from p to x
//         Vector3f p_deviation = (Vector3f::dot(ray.getDirection(), N) < 0) ?
//                 p + N * EPSILON :
//                 p - N * EPSILON ;

//         Ray ray_pTox(p_deviation, ws);
//         //If the ray is not blocked in the middleff
//         Intersection blocked_point_inter = Scene::intersect(ray_pTox);
//         if (abs(distance_pTox - blocked_point_inter.distance < 0.01 ))
//         {
//             L_dir = emit * shade_point_inter.m->eval(wo, ws, N) * Vector3f::dot(ws, N) * Vector3f::dot(-ws, NN) / (distance_pTox * distance_pTox * pdf_light);
//         }
//         //Test Russian Roulette with probability RussianRouolette
//         float ksi = get_random_float();
//         if (ksi < RussianRoulette)
//         {
//             //wi=sample(wo,N)
//             Vector3f wi = (shade_point_inter.m->sample(wo, N)).normalized();
//             //Trace a ray r(p,wi)
//             Ray ray_pTowi(p_deviation, wi);
//             //If ray r hit a non-emitting object at q
//             Intersection bounce_point_inter = Scene::intersect(ray_pTowi);
//             if (bounce_point_inter.isIntersected && !bounce_point_inter.m->hasEmission())
//             {
//                 float pdf = shade_point_inter.m->pdf(wo, wi, N);
// 				if(pdf> EPSILON)
// 					L_indir = castRay(ray_pTowi, depth + 1) * shade_point_inter.m->eval(wo, wi, N) * Vector3f::dot(wi, N) / (pdf *RussianRoulette);
//             }
//         }
//         hitColor = shade_point_inter.m->getEmission() + L_dir + L_indir;
//     }
//     return hitColor;
// }


// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray& ray, int depth) const
{
    //创建变量以储存直接和间接光照计算值
    Vector3f dir = { 0.0,0.0,0.0 };
    Vector3f indir = { 0.0,0.0,0.0 };
    //1.判断是否有交点：光线与场景中物体相交？
    Intersection inter = Scene::intersect(ray);
    //如果没交点
    if (!inter.isIntersected) {
        return dir;//return 0,0,0
    }
    //2.ray打到光源了：说明渲染方程只用算前面的自发光项，因此直接返回材质的自发光项
    if (inter.m->hasEmission()) {
        if (depth == 0) {//第一次打到光
            return inter.m->getEmission();
        }
        else return dir;//弹射打到光，直接返回0，0.0
 
    }
    //3.ray打到物体：这个时候才开始进行伪代码后面的步骤
    
    //对场景中的光源进行采样，得到采样点light_pos和pdf_light
    Intersection light_pos;
    float pdf_light = 0.0f;
    sampleLight(light_pos, pdf_light);
   
    //3.1计算直接光照
 
    //物体的一些参数
    Vector3f p = inter.coords;
    Vector3f N = inter.normal.normalized();
    Vector3f wo = ray.getDirection();//物体指向场景
    //光源的一些参数
    Vector3f xx = light_pos.coords;
    Vector3f NN = light_pos.normal.normalized();
    Vector3f ws = (p - xx).normalized();//光源指向物体
    float dis = (p - xx).length();//二者距离
    float dis2 = Vector3f::dot((p - xx), (p - xx));
 
    //判断光源与物体间是否有遮挡：
    //发出一条射线，方向为ws 光源xx -> 物体p
    Ray light_to_obj(xx, ws);//Ray(orig,dir)
    Intersection light_to_scene = Scene::intersect(light_to_obj);
    //假如dis>light_to_scene.distance就说明有遮挡，那么反着给条件即可：
    if (light_to_scene.isIntersected && (light_to_scene.distance - dis > -EPSILON)) {//没有遮挡
        //为了更贴近伪代码，先设定一些参数
        Vector3f L_i = light_pos.emit;//光强
        Vector3f f_r = inter.m->eval(wo, -ws, N);//材质，课上说了，BRDF==材质，ws不参与计算
        float cos_theta = Vector3f::dot(-ws, N);//物体夹角
        float cos_theta_l = Vector3f::dot(ws, NN);//光源夹角
        dir = L_i * f_r * cos_theta * cos_theta_l / dis2 / pdf_light;
    }
 
    //3.2间接光照
    
    //俄罗斯轮盘赌
    //Scene.hpp中已经定义了P_RR:RussianRoulette=0.8
    float ksi = get_random_float();//随机取[0,1]
    if (ksi < RussianRoulette) {
        //计算间接光照
        
        //随机生成一个wi方向
        Vector3f wi = inter.m->sample(wo, N).normalized();//这里的wi其实没参与计算，返回的是一个随机的方向
        Ray r(p, wi);
        Intersection obj_to_scene = Scene::intersect(r);
        //击中了物体&&物体不是光源
        if (obj_to_scene.isIntersected && !obj_to_scene.m->hasEmission()) {
            Vector3f f_r = inter.m->eval(wo, wi, N);//wo不参与计算
            float cos_theta = Vector3f::dot(wi, N);
            float pdf_hemi = inter.m->pdf(wo, wi, N);
            indir = castRay(r, depth + 1) * f_r * cos_theta / pdf_hemi / RussianRoulette;
        }
    }
    return dir + indir;
}