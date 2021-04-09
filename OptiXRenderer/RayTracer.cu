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
rtDeclareVariable(intersectionData, intersectData, attribute intersectData, );
//rtDeclareVariable(ShadowPayload, shadowPayload, rtPayload, );
rtDeclareVariable(float1, t, rtIntersectionDistance, );
rtDeclareVariable(float3, eye, , );
// depth from Renderer.cpp
rtDeclareVariable(int, depth, , );

RT_PROGRAM void closestHit()
{
    // TODO: calculate the color using the Blinn-Phong reflection model
    float epsilon = .001f;
    float3 result = attrib.ambient + attrib.emission;
    //printf("%f, %f, %f\n", attrib.ambient.x, attrib.ambient.y, attrib.ambient.z);
    //float3 result = attrib.ambient;
     
    // QUESTION: is the HALF-ANGLE: H = normalize(L + V), where L is the direction from hitpoint to light, and V is dir from hitPt to eye?
    for (int i = 0; i < plights.size(); i++) {
        //cast shadow ray
        float3 lightVec = (plights[i].light_pos - intersectData.hitPoint);
        float distToLight = length(lightVec);
        lightVec = normalize(lightVec);

        // create shadow ray and cast it
    	Ray shadowRay = make_Ray(intersectData.hitPoint, lightVec, 1, epsilon, distToLight);
        ShadowPayload shadowPayload; 
        shadowPayload.isVisible = true;
    	rtTrace(root, shadowRay, shadowPayload); 

        //rtPrintf("%d", shadowPayload.isVisible);

        float3 half_angle = normalize(lightVec + normalize(intersectData.rayOrig - intersectData.hitPoint));
        //float3 half_angle = normalize((-intersectData.rayDir) + (eye - intersectData.hitPoint));

        if (shadowPayload.isVisible) {
            //result = (half_angle);
            result += (plights[i].light_color /
                (plights[i].attenuation.constant + plights[i].attenuation.linear * distToLight +
                    plights[i].attenuation.quadratic * powf(distToLight, 2.0f))) * 
                (attrib.diffuse * fmaxf(dot(intersectData.hitPointNormal,
                        normalize(lightVec)), .0f) +
                        attrib.specular * powf(fmaxf(dot( intersectData.hitPointNormal, half_angle), .0f),
                            attrib.shininess));
        }
    }

    for (int i = 0; i < dlights.size(); i++) {
        //cast shadow ray
    	float distToLight = RT_DEFAULT_MAX;
    	Ray shadowRay = make_Ray(intersectData.hitPoint, -dlights[i].light_dir, 1, epsilon, distToLight);
        ShadowPayload shadowPayload; 
        shadowPayload.isVisible = true;
    	rtTrace(root, shadowRay, shadowPayload); 

        float3 half_angle = normalize(normalize(-dlights[i].light_dir) + (intersectData.rayOrig - intersectData.hitPoint));
        //float3 half_angle = normalize(-intersectData.rayDir + (eye - intersectData.hitPoint));
        if (shadowPayload.isVisible) {
            result += (dlights[i].light_color /
                (dlights[i].attenuation.constant + dlights[i].attenuation.linear /** distToLight*/ +
                    dlights[i].attenuation.quadratic /** powf(distToLight, 2.0f)*/)) * (attrib.diffuse * fmaxf(
                        dot(intersectData.hitPointNormal, normalize(-dlights[i].light_dir)), .0f) +
                        attrib.specular * powf(fmaxf(dot(intersectData.hitPointNormal, half_angle), .0f),
                            attrib.shininess));
        }
    }

    // pass the new ray dir and reflection dir into payload 
    // to be used in rayGeneration do-While loop: 
    payload.rayOrigin = intersectData.hitPoint /*+ epsilon * intersectData.hitPointNormal*/;
    payload.rayDir = intersectData.reflectDir;

    if (payload.depth == depth) {
        payload.radiance = result;
        payload.spec = attrib.specular;
    }
    //else payload.radiance = pow(attrib.specular, (depth - (payload.depth + 1)));
    else {
        //float specularExp = depth - (payload.depth + 1); 
        //payload.radiance = make_float3(powf(attrib.specular.x, specularExp), 
        //    powf(attrib.specular.y, specularExp),
        //    powf(attrib.specular.z, specularExp))  * result;
        payload.radiance = result;
        payload.spec *= attrib.specular;
    }
    --payload.depth;
}