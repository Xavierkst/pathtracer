#pragma once

#include <optixu/optixu_math_namespace.h>
#include "Geometries.h"

/**
 * Structures describing different payloads should be defined here.
 */
// The data that describes a ray's journey as it "walks" its path
struct Payload
{
    optix::float3 radiance, throughput, origin, dir;
    unsigned int depth, seed;
    bool done;
};


struct ShadowPayload
{
    int isVisible;
    objectType objType; 
    float3 intersectPt;
};