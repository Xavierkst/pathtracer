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
rtDeclareVariable(ShadowPayload, shadowPayload, rtPayload, );
rtDeclareVariable(float1, t, rtIntersectionDistance, );
rtDeclareVariable(float3, eye, , );

RT_PROGRAM void closestHit()
{
    // TODO: calculate the color using the Blinn-Phong reflection model

    //float3 result = make_float3(0, 1, 0);
    float3 result = attrib.ambient + attrib.emission;
     
    // QUESTION: is the HALF-ANGLE: H = normalize(L + V), where L is the direction from hitpoint to light, and V is dir from hitPt to eye?
    //for (int i = 0; i < plights.size(); i++) {
    //    float3 half_angle = normalize((plights[i].light_pos - payload.hitPoint) + (eye - payload.hitPoint));

    //    result +=  plights[i].light_color / 
    //        (plights[i].attenuation.constant + plights[i].attenuation.linear + 
    //            plights[i].attenuation.quadratic) * (attrib.diffuse * fmaxf(dot(payload.hitPointNormal, (plights[i].light_pos - payload.hitPoint)), .0f) + 
    //                attrib.specular * pow(fmaxf(dot(payload.hitPointNormal, half_angle), .0f), attrib.shininess)); ;
    //}
    //for (int i = 0; i < dlights.size(); i++) {
    //    float3 half_angle = normalize(-dlights[i].light_dir + (eye - payload.hitPoint));

    //    result += dlights[i].light_color /
    //        (dlights[i].attenuation.constant + dlights[i].attenuation.linear +
    //            dlights[i].attenuation.quadratic) * (attrib.diffuse * fmaxf(dot(payload.hitPointNormal, -dlights[i].light_dir), .0f) + 
    //                attrib.specular * pow(fmaxf(dot(payload.hitPointNormal, half_angle), .0f), attrib.shininess));
    //    //result += 
    //}
    //result += shadowPayload.isVisible * result;
    payload.radiance = result;
}