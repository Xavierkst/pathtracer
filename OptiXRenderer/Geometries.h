#pragma once

#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>

/**
 * Structures describing different geometries should be defined here.
 */
struct Attributes
{
    // TODO: define the attributes structure
    optix::float3 ambient;
    float shininess; 
    optix::float3 diffuse;
    optix::float3 specular;
    optix::float3 emission;
};

struct Triangle
{
    // TODO: define the triangle structure
    int maxverts; 
    int maxvertnorms;
    optix::float3 vertex;
    optix::float3 normal;
    int triIndices[3];
    int normIndices[3];

    // include triangle materials here: 
    Attributes attributes;

    Triangle() {
        vertex = optix::make_float3(.0f);
        normal = optix::make_float3(.0f);
        attributes.ambient = optix::make_float3(.0f);
    }
};

struct Sphere
{
    // TODO: define the sphere structure
    optix::float3 position;
    float radius;

    Attributes attributes;
    // include Sphere materials here:

    Sphere() {
        position = optix::make_float3(.0f);
        attributes.ambient = optix::make_float3(.0f);
    }
};

