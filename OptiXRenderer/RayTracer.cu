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
rtDeclareVariable(float1, t, rtIntersectionDistance, );
rtDeclareVariable(int, depth, , ); // depth from Renderer.cpp

RT_PROGRAM void closestHit()
{
    // TODO: calculate the color using the Blinn-Phong reflection model
    float epsilon = .001f;
    float3 result = attrib.ambient + attrib.emission;
     
    for (int i = 0; i < plights.size(); i++) {
        //cast shadow ray
        float3 lightVec = (plights[i].light_pos - intersectData.hitPoint);
        float distToLight = length(lightVec);
        lightVec = normalize(lightVec);

        // create shadow ray and cast it
        float3 shadowRayOrigin = intersectData.hitPoint + intersectData.hitPointNormal * epsilon;
    	Ray shadowRay = make_Ray(shadowRayOrigin, lightVec, 1, epsilon, distToLight);
        ShadowPayload shadowPayload; 
        shadowPayload.isVisible = true;
    	rtTrace(root, shadowRay, shadowPayload); 

        // half angle = normalize of L + V, L (hitPt to light), 
        // V (view dir, eye / prev hitPt to hitPt)
        float3 half_angle = normalize(lightVec + normalize(intersectData.rayOrig - intersectData.hitPoint));

        // compute blinn-phong
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
        float3 shadowRayOrigin = intersectData.hitPoint + intersectData.hitPointNormal * epsilon;
    	Ray shadowRay = make_Ray(shadowRayOrigin, normalize(dlights[i].light_dir), 1, epsilon, distToLight);
    	//Ray shadowRay = make_Ray(intersectData.hitPoint, -dlights[i].light_dir, 1, epsilon, distToLight);
        ShadowPayload shadowPayload; 
        shadowPayload.isVisible = true;
    	rtTrace(root, shadowRay, shadowPayload); 

        // half angle = normalize of L + V, L (hitPt to light), 
        float3 half_angle = normalize(normalize(dlights[i].light_dir) + (intersectData.rayOrig - intersectData.hitPoint));

        // compute blinn-phong
        if (shadowPayload.isVisible) {
            result += (dlights[i].light_color * (attrib.diffuse * fmaxf(
                        dot(intersectData.hitPointNormal, normalize(dlights[i].light_dir)), .0f) +
                        attrib.specular * powf(fmaxf(dot(intersectData.hitPointNormal, half_angle), .0f),
                            attrib.shininess)));
        }
    }
    // pass the new ray dir and reflection dir into payload 
    // to be used in rayGeneration do-While loop: 
    payload.rayOrigin = intersectData.hitPoint + intersectData.hitPointNormal * epsilon;
    payload.rayDir = intersectData.reflectDir;

    if (payload.depth == depth) {
        payload.radiance = result;
        payload.spec = attrib.specular;
    }
    else {
        payload.radiance = payload.spec * result;
        // accumulate specular: 
        // r_1 + S_1 * r2 + S_1 * S_2 * r3 + ...
        payload.spec *= attrib.specular;
    }
}