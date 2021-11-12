#pragma once

#include <optixu/optixu_math_namespace.h>

enum sampling_type { HEMISPHERE_SAMPLING, COSINE_SAMPLING, BRDF_SAMPLING};

enum nextEventEst { OFF, ON, MIS };

struct Config
{
    // Camera 
    optix::float3 w, u, v, eye; // w, u, v: orthonormal basis of camera; eye: eye location 
    optix::float2 hSize, tanHFov; // hSize: half size; tanHFov: tan of 0.5 * fov

    float gamma; 
    // Ray tracing 
    unsigned int maxDepth;
    unsigned int next_event_est;
    float epsilon;
    unsigned int russian_roul;
    
    Config()
    {
        gamma = 1.0f;
        w = optix::make_float3(1.f, 0.f, 0.f);
        u = optix::make_float3(0.f, 1.f, 0.f);
        v = optix::make_float3(0.f, 0.f, 1.f);
        eye = optix::make_float3(0.f);
        hSize = optix::make_float2(100.f); // half the width and height of the resolution
        tanHFov = optix::make_float2(tanf(0.25f * M_PIf));

        maxDepth = 5;
        next_event_est = 0;
        russian_roul = 0;
        epsilon = 0.0001f;
    }
};