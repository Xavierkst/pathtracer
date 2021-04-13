#pragma once

#include <optixu/optixu_math_namespace.h>
#include "Geometries.h"

/**
 * Structures describing different payloads should be defined here.
 */

struct Payload
{
    optix::float3 radiance;
    optix::float3 spec;
    bool done;
    // TODO: add more variable to payload if you need to
    int depth; // recursion depth
    optix::float3 rayOrigin;
    optix::float3 rayDir;
};

struct ShadowPayload
{
    int isVisible;
};