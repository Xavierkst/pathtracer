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
    for (int i = 0; i < qlights.size(); i++)
    {

	float3 points[] = {qlights[i].tri1.v1,qlights[i].tri1.v2,
			   qlights[i].tri1.v3,qlights[i].tri2.v2};
	float3 phi = make_float3(0);
	float3 k_d = mv.diffuse;
	float3 L_i = qlights[i].color;
	float3 normal = attrib.normal;
	float3 r = attrib.intersection;

	for (int j = 0; j < 4; j++) {
	    float3 p1 = points[j];
	    float3 p2 = points[(j + 1)%4];

	    float theta = acos(dot(normalize(p1 - r), normalize(p2 - r)));
	    float3 gamma = normalize(cross(p1 - r, p2 - r));

	    phi += theta * gamma;
	}
	phi *= .5f;

	result += ((k_d/M_PIf)*L_i)*(dot(phi, normal));


        // Shoot a shadow to determin whether the object is in shadow
	/*
        ShadowPayload shadowPayload;
        shadowPayload.isVisible = true;
        Ray shadowRay = make_Ray(attrib.intersection + lightDir * cf.epsilon, 
            lightDir, 1, cf.epsilon, lightDist);
        rtTrace(root, shadowRay, shadowPayload);

        // If not in shadow
        if (shadowPayload.isVisible)
        {
        }
	*/
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
