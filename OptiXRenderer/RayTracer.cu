#include <optix.h>
#include <optix_device.h>
#include <optixu/optixu_math_namespace.h>
#include "random.h"

#include "Payloads.h"
#include "Geometries.h"
#include "Light.h"

using namespace optix;

// Declare light buffers
rtBuffer<PointLight> plights;
rtBuffer<DirectionalLight> dlights;

// Declare variables
rtDeclareVariable(Payload, payload, rtPayload, );
rtDeclareVariable(rtObject, root, , );

// Declare attibutes 
rtDeclareVariable(Attributes, attrib, attribute attrib, );

RT_PROGRAM void closestHit()
{
    // TODO: calculate the color using the Blinn-Phong reflection model

    //float3 result = make_float3(0, 1, 0);
    float3 result = attrib.ambient;
    result += dlights[0].light_color;
    for (int i = 0; i < plights.size(); i++) {
        result += plights[i].light_color;
    }
    payload.radiance = result;
}