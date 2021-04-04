#pragma once

#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>

#include "Geometries.h"
#include "Light.h"
#include "Camera.h"

struct Scene
{
    // Info about the output image
    std::string outputFilename;
    unsigned int width, height;

    std::string integratorName;

    std::vector<optix::float3> vertices;
    std::vector<optix::float3> normals;
    // a vector of vertex + normal pairs
    // first == vertex, second == normal
    std::vector<std::pair<optix::float3, optix::float3>> vertexNormals;

    std::vector<Triangle> triangles;
    std::vector<Sphere> spheres;
    unsigned int maxverts;
    unsigned int maxvertnorms;

    std::vector<DirectionalLight> dlights;
    std::vector<PointLight> plights;

    // TODO: add other variables that you need here
    Camera camera;
    unsigned int depth;

    Scene()
    {
        depth = 5;
        outputFilename = "raytrace.png";
        integratorName = "raytracer";
    }
};