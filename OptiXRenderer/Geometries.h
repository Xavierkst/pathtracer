#pragma once

#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>

/**
 * Structures describing different geometries should be defined here.
 */

enum objectType { OBJECT, LIGHT };

enum brdfType { MOD_PHONG, GGX/*, VOLUMETRIC */};

enum materialType { DIFFUSE, GLASS, VOLUMETRIC };

struct MaterialValue
{
    optix::float3 ambient, diffuse, specular, emission;
    float shininess, roughness;
    brdfType brdf_type;    
    float ior;
    materialType matType;
    float sigma_a, sigma_s, sigma_t, g;
};

struct Triangle
{
    optix::float3 v1, v2, v3, normal; // transformed
    MaterialValue mv;
    objectType objType;
};

struct Sphere
{
    optix::Matrix4x4 trans;
    MaterialValue mv;
    objectType objType;
};

struct Attributes
{
    optix::float3 intersection, normal, wo, prev_intersection;
    MaterialValue mv;
    objectType objType;
};
