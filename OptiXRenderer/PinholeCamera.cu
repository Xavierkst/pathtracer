#include <optix.h>
#include <optix_device.h>
#include <optixu/optixu_math_namespace.h>

#include "Payloads.h"
#include "Camera.h"

using namespace optix;

rtBuffer<float3, 2> resultBuffer; // used to store the render result

rtDeclareVariable(rtObject, root, , ); // Optix graph

rtDeclareVariable(uint2, launchIndex, rtLaunchIndex, ); // a 2d index (x, y)

rtDeclareVariable(int1, frameID, , );

// Camera info 

// TODO:: delcare camera variables here
rtDeclareVariable(float3, eye, , );
rtDeclareVariable(float3, U, , );
rtDeclareVariable(float3, V, , );
rtDeclareVariable(float3, W, , );
rtDeclareVariable(float, fovy, , );
rtDeclareVariable(int, width, , );
rtDeclareVariable(int, height, , );
rtDeclareVariable(int, depth, , );

//rtPrintf("%d", resultBuffer.size());

RT_PROGRAM void generateRays()
{
    float3 result = make_float3(0.f);
    // TODO: calculate the ray direction (change the following lines)
    float3 origin = eye;  // origin should be pos of camera
    float aspectRatio = (float) width / (float)height;
    float alpha = ((2.0f * ((float)launchIndex.x + 0.5f) / (float)width) - 1.0f) * tan(fovy / 2.0f) * aspectRatio;
    float beta = ((2.0f * ((float)launchIndex.y + 0.5f) / (float)height) - 1.0f) * tan(fovy / 2.0f);

    float3 dir = normalize(alpha * U + beta * V - W);

    float epsilon = 0.001f; 
    Payload payload;
    payload.depth = depth;
    // TODO: modify the following lines if you need
    // Shoot a ray to compute the color of the current pixel
    //Ray ray = make_Ray(origin, dir, 0, epsilon, RT_DEFAULT_MAX);
    //rtTrace(root, ray, payload);
    //result = payload.radiance;


    do {
        Ray ray = make_Ray(origin, dir, 0, epsilon, RT_DEFAULT_MAX);
        rtTrace(root, ray, payload);
        result += payload.radiance;
        // set up for next ray cast
        origin = payload.hitPoint; 
        dir = payload.dir;
        //--payload.depth;
    } while (!payload.done && payload.depth > 0);

    // Write the result
    resultBuffer[launchIndex] = result;
}