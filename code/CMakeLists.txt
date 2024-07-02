CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(PathTracing CXX)

IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF()

ADD_SUBDIRECTORY(deps/vecmath)

SET(PA1_SOURCES
        src/main.cpp
        src/BVH.cpp
        src/Renderer.cpp
        src/Scene.cpp)

SET(PA1_INCLUDES
        include/Bounds3.hpp
        include/BVH.hpp
        include/global.hpp
        include/Intersection.hpp
        include/Light.hpp
        include/Material.hpp
        include/OBJ_Loader.hpp
        include/object3d.hpp
        include/plane.hpp
        include/ray.hpp
        include/Renderer.hpp
        include/Scene.hpp
        include/sphere.hpp
        include/triangle.hpp
        )

SET(CMAKE_CXX_STANDARD 17)

ADD_EXECUTABLE(${PROJECT_NAME} ${PA1_SOURCES} ${PA1_INCLUDES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} vecmath)
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE include)
