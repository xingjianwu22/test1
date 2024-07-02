#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include <mutex>
#include <thread>

inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00016;
std::mutex mtx;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    // change the spp value to change sample ammount
    int spp = 256;                 // spp指每个pixel会采样的次数
    std::cout << "SPP: " << spp << "\n";
    //int m = 0;
    // for (uint32_t j = 0; j < scene.height; ++j) {
        // for (uint32_t i = 0; i < scene.width; ++i) {
        //     for (int k = 0; k < spp; k++) {
        //         // 超采样抗锯齿（SSAA） 随机超采样（Random Sampling SSAA）
        //         // 生成随机偏差 [-0.5, 0.5]
        //         float offset_x = ((rand() / (float)RAND_MAX) - 0.5f);
        //         float offset_y = ((rand() / (float)RAND_MAX) - 0.5f);

        //         // 计算子像素位置
        //         float x = (2 * (i + 0.5f + offset_x) / (float)scene.width - 1) * imageAspectRatio * scale;
        //         float y = (1 - 2 * (j + 0.5f + offset_y) / (float)scene.height) * scale;

        //         // 计算光线方向
        //         Vector3f dir = Vector3f(-x, y, 1).normalized();

        //         // 发射光线并累加颜色值
        //         framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
        //     }
        //     m++;
    //     }
    //     UpdateProgress(j / (float)scene.height);
    // }
    // UpdateProgress(1.f);

    int process = 0;
    const int thred = 20;
    int times = scene.height / thred;  // 每个线程的行数
    std::thread th[thred];

    auto castRayMultiThread = [&](uint32_t y_min, uint32_t y_max)
    {
        for(uint32_t j = y_min; j < y_max; j++)
        {
            
            int m = j * scene.width;
            for (uint32_t i = 0; i < scene.width; ++i) {
                for (int k = 0; k < spp; k++) {
                    // 超采样抗锯齿（SSAA） 随机超采样（Random Sampling SSAA）
                    // 生成随机偏差 [-0.5, 0.5]
                    float offset_x = ((rand() / (float)RAND_MAX) - 0.5f);
                    float offset_y = ((rand() / (float)RAND_MAX) - 0.5f);

                    // 计算子像素位置
                    float x = (2 * (i + 0.5f + offset_x) / (float)scene.width - 1) * imageAspectRatio * scale;
                    float y = (1 - 2 * (j + 0.5f + offset_y) / (float)scene.height) * scale;

                    // 计算光线方向
                    Vector3f dir = Vector3f(-x, y, 1).normalized();

                    // 发射光线并累加颜色值
                    framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;
                }
                m++;
            }
            mtx.lock();
            process++;
            UpdateProgress(1.0 * process / scene.height);
            mtx.unlock();
        }
    };


    //分行进行路径追踪
    for (int i = 0; i < thred; i++) {//从第0行出发，一共有0~by-1行
        th[i] = std::thread(castRayMultiThread, i * times, (i + 1) * times);
    }
    //每个线程执行join
    for (int i = 0; i < thred; i++) {
        th[i].join();
    }
    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        // 伽马校正
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x()), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y()), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z()), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
