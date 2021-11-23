#include <optix.h>
#include <optix_device.h>
#include <optixu/optixu_math_namespace.h>
#include "random.h"

#include "Payloads.h"
#include "Config.h"
#include "Light.h"

using namespace optix;

rtBuffer<float3, 2> resultBuffer; // used to store the render result
rtDeclareVariable(rtObject, root, , ); // Optix graph
rtDeclareVariable(uint2, launchIndex, rtLaunchIndex, ); // a 2d index (x, y)
rtDeclareVariable(int1, frameID, , );
rtBuffer<Config> config; // Config
rtDeclareVariable(int, samples_per_pixel, , );
rtDeclareVariable(uint, next_event_est, , );

RT_PROGRAM void generateRays()
{
    size_t2 resultSize = resultBuffer.size();
    unsigned int index = launchIndex.x * resultSize.y + launchIndex.y;
    unsigned int seed = tea<16>(index * frameID.x, 0);
    Config cf = config[0];
    float3 result = make_float3(0.f);
    float2 xy = make_float2(launchIndex);
    // xy.x += frameID.x == 1 ? 0.5f : rnd(seed);
    // xy.y += frameID.x == 1 ? 0.5f : rnd(seed);
    float2 ab;
    float3 dir;
    float3 origin;
    Payload payload;
    int i = 0;
    // cf.maxDepth = RT_DEFAULT_MAX;
    // Cast samples_per_pixel number of rays thru pixel xy
    for (int j = 0; j < samples_per_pixel; ++j) {
        // Prepare new payload for each sample
        payload.radiance = make_float3(.0f);
        payload.throughput = make_float3(1.0f);
        payload.depth = 0;
        payload.done = false;
        // Compute the ray direction: 
        xy = make_float2(launchIndex);
        // for th very first sample, we keep it at the center of the pixel
        // For every subsequent sample, jitter rays entering pixel
        xy.x += (j == 0) ? 0.5f : rnd(seed);
        xy.y += (j == 0) ? 0.5f : rnd(seed);
        ab = cf.tanHFov * (xy - cf.hSize) / cf.hSize; // calculates NDC coordinates -1 to +1 
        origin = cf.eye;
        dir = normalize(ab.x * cf.u + ab.y * cf.v - cf.w); // ray direction

        // For each pixel sample, we trace the path up to a depth D == cf.maxDepth
        do
        {
            payload.seed = tea<16>(index * frameID.x, i++);
            // Cast primary ray into the scene
            Ray ray = make_Ray(origin, dir, 0, cf.epsilon, RT_DEFAULT_MAX);
            rtTrace(root, ray, payload); // Goes to intersection program for Sphere and Tri
            // Accumulate radiance
            result += payload.radiance;
            payload.radiance = make_float3(0.f);
            // Continue "walking" thru the scene from one hit point to another
            origin = payload.origin; 
            dir = payload.dir;
        } while (!payload.done && payload.depth != cf.maxDepth);
    }
    
    // average out the results 
    result = (result / samples_per_pixel);
    result = make_float3(powf(result.x, 1.0f / cf.gamma), powf(result.y, 1.0f / cf.gamma), powf(result.z, 1.0f / cf.gamma));
    
    if (frameID.x == 1) 
        resultBuffer[launchIndex] = result;
    else
    {
        float u = 1.0f / (float)frameID.x;
        float3 oldResult = resultBuffer[launchIndex];
        resultBuffer[launchIndex] = lerp(oldResult, result, u);
    }
}