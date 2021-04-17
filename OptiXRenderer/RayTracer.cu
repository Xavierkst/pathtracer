#include <optix.h>
#include <optix_device.h>
#include <optixu/optixu_math_namespace.h>
#include "random.h"

#include "Payloads.h"
#include "Geometries.h"
#include "Light.h"
#include "Config.h"

using namespace optix;

// Declare light buffers
rtBuffer<PointLight> plights;
rtBuffer<DirectionalLight> dlights;
rtBuffer<QuadLight> qlights;

// Declare variables
rtDeclareVariable(Payload, payload, rtPayload, );
rtDeclareVariable(rtObject, root, , );

rtBuffer<Config> config; // Config

// Declare attibutes 
rtDeclareVariable(Attributes, attrib, attribute attrib, );

RT_PROGRAM void closestHit()
{
    MaterialValue mv = attrib.mv;
    Config cf = config[0];

    float3 result = mv.ambient + mv.emission;

    // Calculate the direct illumination of point lights
    for (int i = 0; i < plights.size(); i++)
    {
        // Shoot a shadow to determin whether the object is in shadow
        float3 lightDir = normalize(plights[i].location - attrib.intersection);
        float lightDist = length(plights[i].location - attrib.intersection);
        ShadowPayload shadowPayload;
        shadowPayload.isVisible = true;
        Ray shadowRay = make_Ray(attrib.intersection + lightDir * cf.epsilon, 
            lightDir, 1, cf.epsilon, lightDist);
        rtTrace(root, shadowRay, shadowPayload);

        // If not in shadow
        if (shadowPayload.isVisible)
        {
            float3 H = normalize(lightDir + attrib.wo);
            float att = dot(plights[i].attenuation, make_float3(1, lightDist, lightDist * lightDist));
            float3 I = mv.diffuse * fmaxf(dot(attrib.normal, lightDir), 0);
            I += mv.specular * pow(fmaxf(dot(attrib.normal, H), 0), mv.shininess);
            I *= plights[i].color / att;
            result += I;
        }
    }

    // Calculate the direct illumination of directional lights
    for (int i = 0; i < dlights.size(); i++)
    {
        // Shoot a shadow to determin whether the object is in shadow
        float3 lightDir = dlights[i].direction;
        float lightDist = RT_DEFAULT_MAX;
        ShadowPayload shadowPayload;
        shadowPayload.isVisible = true;
        Ray shadowRay = make_Ray(attrib.intersection + lightDir * cf.epsilon, 
            lightDir, 1, cf.epsilon, lightDist);
        rtTrace(root, shadowRay, shadowPayload);

        // If not in shadow
        if (shadowPayload.isVisible)
        {
            float3 H = normalize(lightDir + attrib.wo);
            float3 I = mv.diffuse * fmaxf(dot(attrib.normal, lightDir), 0);
            I += mv.specular * pow(fmaxf(dot(attrib.normal, H), 0), mv.shininess);
            I *= dlights[i].color;
            result += I;
        }
    }

    // Another for-loop here to calculate contribution of quadLights
    // eg. for (int i = 0; i < qlights.size(); ++i) {}
    for (int i = 0; i < qlights.size(); ++i) {
        float3 f_brdf = mv.diffuse / M_PIf;// brdf function 
        float3 hitPt = attrib.intersection;
        float3 hitPtNormal = attrib.normal;

        float3 a = qlights[i].tri1->v1;
        float3 b = qlights[i].tri1->v2;
        float3 c = qlights[i].tri2->v2;
        float3 d = qlights[i].tri1->v3;


        float theta_1 = acosf(dot(normalize(a - hitPt), normalize(b - hitPt)));
        float theta_2 = acosf(dot(normalize(b - hitPt), normalize(c - hitPt)));
        float theta_3 = acosf(dot(normalize(c - hitPt), normalize(d - hitPt)));
        float theta_4 = acosf(dot(normalize(d - hitPt), normalize(a - hitPt)));

        float3 gamma_1 = normalize(cross((a - hitPt), (b - hitPt)));
        float3 gamma_2 = normalize(cross((b - hitPt), (c - hitPt)));
        float3 gamma_3 = normalize(cross((c - hitPt), (d - hitPt)));
        float3 gamma_4 = normalize(cross((d - hitPt), (a - hitPt)));

        float3 irradiance_vec = 0.5f * (theta_1 * gamma_1 + 
            theta_2 * gamma_2 + theta_3 * gamma_3 * theta_4 * gamma_4);

        float3 dir_radiance = f_brdf * qlights[i].color * dot(irradiance_vec, hitPtNormal);
        result += dir_radiance;
    }

    // Compute the final radiance
    payload.radiance = result * payload.throughput;

    // Calculate reflection
    if (length(mv.specular) > 0)
    {
        // Set origin and dir for tracing the reflection ray
        payload.origin = attrib.intersection;
        payload.dir = reflect(-attrib.wo, attrib.normal); // mirror reflection

        payload.depth++;
        payload.throughput *= mv.specular;
    }
    else
    {
        payload.done = true;
    }
}