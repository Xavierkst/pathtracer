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

rtDeclareVariable(uint, light_samples, , );
rtDeclareVariable(uint, light_stratify, , );

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
     //eg. for (int i = 0; i < qlights.size(); ++i) {}
    // if no light samples, do analytical direct:
    if (light_samples == 0) {
        for (int i = 0; i < qlights.size(); ++i) {
            float3 f_brdf = mv.diffuse / M_PIf;// brdf function 
            float3 hitPt = attrib.intersection;
            float3 hitPtNormal = attrib.normal;

            float3 a = qlights[i].tri1.v1;
            float3 b = qlights[i].tri1.v2;
            float3 c = qlights[i].tri2.v2;
            float3 d = qlights[i].tri1.v3;

            float3 points[] = { qlights[i].tri1.v1, qlights[i].tri1.v2, qlights[i].tri2.v2, qlights[i].tri1.v3 };

            float3 p1 = points[0]; float3 p2 = points[1]; float3 p3 = points[2]; float3 p4 = points[3];
            float theta_1 = acosf(dot(normalize(p1 - hitPt), normalize(p2 - hitPt)));
            float theta_2 = acosf(dot(normalize(p2 - hitPt), normalize(p3 - hitPt)));
            float theta_3 = acosf(dot(normalize(p3 - hitPt), normalize(p4 - hitPt)));
            float theta_4 = acosf(dot(normalize(p4 - hitPt), normalize(p1 - hitPt)));

            float3 gamma_1 = normalize(cross((p1 - hitPt), (p2 - hitPt)));
            float3 gamma_2 = normalize(cross((p2 - hitPt), (p3 - hitPt)));
            float3 gamma_3 = normalize(cross((p3 - hitPt), (p4 - hitPt)));
            float3 gamma_4 = normalize(cross((p4 - hitPt), (p1 - hitPt)));

            float3 irradiance_vec = 0.5f * ((theta_1 * gamma_1) +
                (theta_2 * gamma_2) + (theta_3 * gamma_3) + (theta_4 * gamma_4));

            float3 dir_radiance = f_brdf * qlights[i].color * dot(irradiance_vec, hitPtNormal);
            result += dir_radiance;
        }
    }
    else { // light_samples > 0
        for (int k = 0; k < qlights.size(); ++k) {
            float3 sampled_result = make_float3(.0f);
            // Compute direct lighting equation for w_i_k ray, for k = 1 to N*N
            float3 a = qlights[k].tri1.v1;
            float3 b = qlights[k].tri1.v2;
            float3 c = qlights[k].tri2.v2;
            float3 d = qlights[k].tri1.v3;

            //float3 f_brdf_1 = mv.diffuse / M_PIf;// brdf function 
            //float3 hitPt = attrib.intersection;
            //float3 hitPtNormal = attrib.normal;

            //float3 points[] = { qlights[k].tri1.v1, qlights[k].tri1.v2, qlights[k].tri2.v2, qlights[k].tri1.v3 };

            //float3 p1 = points[0]; float3 p2 = points[1]; float3 p3 = points[2]; float3 p4 = points[3];
            //float theta_1 = acosf(dot(normalize(p1 - hitPt), normalize(p2 - hitPt)));
            //float theta_2 = acosf(dot(normalize(p2 - hitPt), normalize(p3 - hitPt)));
            //float theta_3 = acosf(dot(normalize(p3 - hitPt), normalize(p4 - hitPt)));
            //float theta_4 = acosf(dot(normalize(p4 - hitPt), normalize(p1 - hitPt)));

            //float3 gamma_1 = normalize(cross((p1 - hitPt), (p2 - hitPt)));
            //float3 gamma_2 = normalize(cross((p2 - hitPt), (p3 - hitPt)));
            //float3 gamma_3 = normalize(cross((p3 - hitPt), (p4 - hitPt)));
            //float3 gamma_4 = normalize(cross((p4 - hitPt), (p1 - hitPt)));

            //float3 irradiance_vec = 0.5f * ((theta_1 * gamma_1) + 
            //    (theta_2 * gamma_2) + (theta_3 * gamma_3) + (theta_4 * gamma_4));

            //float3 dir_radiance = f_brdf_1 * qlights[k].color * dot(irradiance_vec, hitPtNormal);
            //sampled_result += dir_radiance;

            float3 ac = c - a;
            float3 ab = b - a;
            float area = length(cross(ab, ac));
            int root_light_samples = (int)sqrtf(light_samples);
            //rtPrintf("root of light samples and stratify: %d %d\n", root_light_samples, light_stratify);
            // check if stratify or random sampling
                // double for loop here 
            for (int i = 0; i < root_light_samples; ++i) {
                for (int j = 0; j < root_light_samples; ++j) {
                    // generate random float vals u1 and u2
                    float u1 = rnd(payload.seed);
                    float u2 = rnd(payload.seed);
                    //rtPrintf("%f %f \n", u1, u2);
                    float3 sampled_light_pos;
                    if (light_stratify) {
                        sampled_light_pos = a + ((j + u1) * (ab / (float)root_light_samples)) +
                            ((i + u2) * (ac / (float)root_light_samples));
                    }
                    else {
                        sampled_light_pos = a + u1 * ab + u2 * ac;
                    }
                    float3 shadow_ray_origin = attrib.intersection /*+ attrib.normal * cf.epsilon*/;
                    float3 shadow_ray_dir = normalize(sampled_light_pos - shadow_ray_origin);
                    float light_dist = length(sampled_light_pos - shadow_ray_origin);
                    Ray shadow_ray = make_Ray(shadow_ray_origin, shadow_ray_dir, 1, cf.epsilon, light_dist - cf.epsilon);

                    ShadowPayload shadow_payload;
                    shadow_payload.isVisible = true;
                    rtTrace(root, shadow_ray, shadow_payload);

                    //rtPrintf("%d", shadow_payload.isVisible);

                    if (shadow_payload.isVisible) {
                        // rendering equation here: 
                        //float3 w_i = sampled_light_pos;
                        float3 f_brdf = (mv.diffuse / M_PIf) +
                            (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                                powf(fmaxf(dot(normalize(reflect(-attrib.wo, attrib.normal)), normalize(sampled_light_pos - shadow_ray_origin)), .0f), mv.shininess));

                        float3 x_prime = sampled_light_pos;
                        float3 x = shadow_ray_origin;
                        float3 n = attrib.normal;
                        //float3 n_light = normalize(qlights[k].tri1.normal);
                        float3 n_light = normalize(cross(ab, ac));
                        //n_light = dot(n_light, normalize(x_prime - x)) > .0f ? n_light : -n_light;

                        float R = length(x - x_prime);

                        // note: normal should point AWAY from the hitpoint, i.e. dot(n_light, x - x_prime) < 0
                        float G = (1.0f / powf(R, 2.0f)) * fmaxf(dot(n, normalize(x_prime - x)), .0f) *
                            (fmaxf(dot(n_light, normalize(x_prime - x)), .0f));

                        sampled_result += f_brdf * G;
                    }

                }
            }
            result += qlights[k].color * sampled_result * (area / (float)light_samples);
        }
    }
    //}

    //result += sampled_result;



    // Compute the final radiance
    payload.radiance = result * payload.throughput;

    // Calculate reflection
    //if (length(mv.specular) > 0)
    //{
    //    // Set origin and dir for tracing the reflection ray
    //    payload.origin = attrib.intersection;
    //    payload.dir = reflect(-attrib.wo, attrib.normal); // mirror reflection

    //    payload.depth++;
    //    payload.throughput *= mv.specular;
    //}
    //else
    //{
    //    payload.done = true;
    //}
    payload.done = true;
}