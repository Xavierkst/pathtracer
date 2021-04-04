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

//rtPrintf("%d", resultBuffer.size());

RT_PROGRAM void generateRays()
{
    size_t2 screen_size = resultBuffer.size();
    //rtPrintf("the total number of pixels is: %d", screen_size);

    float3 result = make_float3(0.f);

    // TODO: calculate the ray direction (change the following lines)
    float3 origin = eye;  // origin should be pos of camera
    
    //float2 d = make_float2(launchIndex) / make_float2(screen_size) * 2.0f - 1.0f;
    float aspectRatio = (float) width / (float)height;
    float alpha = ((2.0f * ((float)launchIndex.x + 0.5f) / (float)width) - 1.0f) * tan(fovy / 2.0f) * aspectRatio;
    float beta = (1.0f - (2.0f * ((float)launchIndex.y + 0.5f) / (float)height) ) * tan(fovy/2.0f);
    //rtPrintf("alpha beta: %f %f", alpha, beta);

    //float3 dir = make_float3(0, 0, 1); // dir should be toward some (i,j)-th cell?
    float3 dir = normalize(alpha * U + beta * V - W);

    float epsilon = 0.001f; 

    // TODO: modify the following lines if you need
    // Shoot a ray to compute the color of the current pixel
    Ray ray = make_Ray(origin, dir, 0, epsilon, RT_DEFAULT_MAX);
    Payload payload;
    rtTrace(root, ray, payload);

    // Write the result
    resultBuffer[launchIndex] = payload.radiance;
}