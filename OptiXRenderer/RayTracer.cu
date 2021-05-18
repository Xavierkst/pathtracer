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
    //payload.done = true;
}

RT_PROGRAM void analyticDirect() {

    MaterialValue mv = attrib.mv;
    Config cf = config[0];

    float3 result = mv.ambient + mv.emission;

    // for-loop here to calculate contribution of quadLights
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

    payload.radiance = result;

    payload.done = true;
}


RT_PROGRAM void direct() {

    MaterialValue mv = attrib.mv;
    Config cf = config[0];

    float3 result = mv.ambient + mv.emission;

    for (int k = 0; k < qlights.size(); ++k) {
        float3 sampled_result = make_float3(0);
        // Compute direct lighting equation for w_i_k ray, for k = 1 to N*N
        float3 a = qlights[k].tri1.v1;
        float3 b = qlights[k].tri1.v2;
        float3 c = qlights[k].tri2.v3;
        float3 d = qlights[k].tri2.v2;

        float3 ac = c - a;
        float3 ab = b - a;
        float area = length(cross(ab, ac));
        int root_light_samples = (int)sqrtf(light_samples);
        // check if stratify or random sampling
        // double for loop here 
        for (int i = 0; i < root_light_samples; ++i) {
            for (int j = 0; j < root_light_samples; ++j) {
                // generate random float vals u1 and u2
                float u1 = rnd(payload.seed);
                float u2 = rnd(payload.seed);

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
    //rtPrintf("throughput val: %f \n", payload.throughput);
    payload.radiance = result;

    payload.done = true;
}

RT_PROGRAM void pathtracer() {

    MaterialValue mv = attrib.mv;
    Config cf = config[0];

    float pdf = 1.0f;
    float3 brdf = make_float3(0);
    float3 result = make_float3(0);
    float3 L_d = make_float3(0);
    float3 L_e = mv.emission;
    float3 r = normalize(reflect(-attrib.wo, attrib.normal));
    float s_mean = (mv.specular.x + mv.specular.y + mv.specular.z)/3.0f;
    float d_mean = (mv.diffuse.x + mv.diffuse.y + mv.diffuse.z)/3.0f;
    float u0 = rnd(payload.seed);
    float u1 = rnd(payload.seed);
    float u2 = rnd(payload.seed);
    float t = 0.0f;

    if (s_mean + d_mean == 0) {
	if (mv.brdf == 1) t = 1.0f;
    }
    else {
	t = (mv.brdf == 0) ? s_mean/(s_mean+d_mean) : fmaxf(0.25f, s_mean/(s_mean+d_mean));
    }

    float theta;
    float phi;

    switch (cf.IS) {
        case 0: 
            phi = 2.0f * M_PIf * u2;
            theta = acosf(u1);
            break;
        case 1: 
            phi = 2.0f * M_PIf * u2;
            theta = acosf(sqrtf(u1));
            break;
        case 2: 
            phi = 2.0f * M_PIf * u2;
            if (u0 > t) 
                theta = acosf(sqrtf(u1)); // theta_diffuse
            else 
                theta = acosf(powf(u1, (1.0f / (mv.shininess + 1.0f)))); // theta_specular
            break;
    }

    float3 sampleVec = make_float3(cosf(phi)*sinf(theta), sinf(phi)*sinf(theta), cosf(theta));

    float3 n = normalize(attrib.normal);
    float3 w = ((u0 <= t) && cf.IS == 2 && mv.brdf == 0) ? r : n;
    float3 a = make_float3(0,1,0);
    a = fabsf(dot(a,w)) > .9f ? make_float3(1,0,0) : a;

    float3 u = normalize(cross(a, w));
    float3 v = cross(w,u);
    if (mv.brdf == 1) {
	float theta_h_sample = atanf((mv.roughness * sqrtf(u1)) / sqrtf(1.0f - u1));
	float phi_h_sample = 2.0f * M_PIf * u2;
	sampleVec = make_float3(cosf(phi_h_sample) * sinf(theta_h_sample), 
	       sinf(phi_h_sample)*sinf(theta_h_sample), 
	       cosf(theta_h_sample));
    }

    float3 wi = make_float3(1.0f); 

    // get randomized new ray dir -- choose a sampling method
    wi = (sampleVec.x * u + sampleVec.y * v + sampleVec.z * w);   

    if (mv.brdf == 1) {
        wi = (reflect(-attrib.wo, wi));
    }

    float3 bruh = make_float3(.0f);

    switch (cf.IS) {
        case 0: 
            brdf = (mv.diffuse / M_PIf) +
                (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                    powf(fmaxf(dot(r, wi), .0f), mv.shininess));
            pdf = 1.0f / (2.0f * M_PIf);
            bruh = (brdf * fmaxf(dot(n, wi), .0f) * (1.0f / pdf));
            break;

        case 1: 
            brdf = (mv.diffuse / M_PIf) +
                (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                    powf(fmaxf(dot(r, wi), .0f), mv.shininess));
            pdf = dot(n, wi) / (M_PIf);
            bruh = (brdf * fmaxf(dot(n, wi), .0f) * (1.0f / pdf));
            break;

        case 2: 
            // check the material whether to use mod-phong or GGX brdf
            if (mv.brdf == 0) {
                brdf = (mv.diffuse / M_PIf) +
		       (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                        powf(fmaxf(dot(r, wi), .0f), mv.shininess));

                pdf = ((1.0f - t) * (fmaxf(dot(n, wi), .0f) / M_PIf)) + 
			t * ((mv.shininess + 1.0f) / (2.0f * M_PIf)) * 
                        powf(fmaxf(dot(r, wi), .0f), mv.shininess);
            }
            else {
                // construct GGX BRDF: 
		float wi_dot_n = dot(wi, n);
		float wo_dot_n = dot(attrib.wo, n);
                if (wi_dot_n > .0f && wo_dot_n > .0f) {
                    float alpha = mv.roughness;
                    float3 h = normalize(wi + attrib.wo); // half angle: 
                    float theta_h = acosf(clamp(dot(h, n),0.0f,1.0f)); // not sure if need to clamp 0
                    // microfacet distribution function, D: 
                    float D = (alpha * alpha) / (M_PIf * powf(cosf(theta_h), 4.0f) *
			       powf((alpha * alpha) + powf(tanf(theta_h), 2.0f), 2.0f));
                    // shadow-masking function, G:  
                    float G_1_wi = (wi_dot_n > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) *
				    powf(tanf(acosf(clamp(dot(wi, n),0.0f,1.0f))), 2.0f))) : .0f;
                    float G_1_wo = (wo_dot_n > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) *
				    powf(tanf(acosf(clamp(dot(attrib.wo, n),0.0f,1.0f))), 2.0f))) : .0f;
                    float G = G_1_wi * G_1_wo;

                    // fresnel function, F:
                    float3 F = mv.specular + (make_float3(1.0f) - mv.specular) * powf(1.0f - dot(wi, h), 5.0f);
                    float3 f_brdf_GGX = (F * G * D) / (4.0f * wi_dot_n * wo_dot_n);
                    brdf = (mv.diffuse / M_PIf) + f_brdf_GGX;

                    pdf = fmaxf(((1 - t) * (wi_dot_n / M_PIf)) + ((t * D * dot(n, h)) / (4.0f * fmaxf(dot(h, wi),0.0f))),0.0f);
                }
                else brdf = make_float3(0.0f); // assume f zero otherwise
            }
            bruh = (brdf * fmaxf(dot(n,wi),0) * (1.0f / pdf));
            break;
    }


    for (int k = 0; k < (cf.NEE ? qlights.size()*cf.NEE : 0); ++k) {
        float3 sampled_result = make_float3(.0f);
	float3 brdf_cum = make_float3(0);
        // Compute direct lighting equation for w_i_k/2 ray, for k = 1 to N*N
	int whichLight = (k >= qlights.size()) ? k%qlights.size() : k;
        float3 a = qlights[whichLight].tri1.v1;
        float3 b = qlights[whichLight].tri1.v2;
        float3 c = qlights[whichLight].tri2.v3;
        float3 d = qlights[whichLight].tri2.v2;

        float3 ac = c - a + cf.epsilon;
        float3 ab = b - a + cf.epsilon;
        float area = length(cross(ab, ac));
        int root_light_samples = (int)sqrtf(light_samples);
        // check if stratify or random sampling
        // double for loop here 
        for (int i = 0; i < root_light_samples; ++i) {
            for (int j = 0; j < root_light_samples; ++j) {

                float3 sampled_light_pos = make_float3(0);
                if (light_stratify) {
                    sampled_light_pos = a + ((j + u1) * (ab / (float)root_light_samples)) +
                        ((i + u2) * (ac / (float)root_light_samples));
                }
                else {
                    sampled_light_pos = a + u1 * ab + u2 * ac;
                }

                float3 shadow_ray_origin = attrib.intersection; //+ attrib.normal * cf.epsilon;
                float3 lightDir = normalize(sampled_light_pos - shadow_ray_origin);
                float light_dist = length(sampled_light_pos - shadow_ray_origin);
                ShadowPayload shadow_payload;
                shadow_payload.isVisible = true;

		if (k >= qlights.size()) {

		    for (int poggers = 0; poggers < 2; poggers++) {

			Triangle tri = poggers ? qlights[whichLight].tri1 : qlights[whichLight].tri2;

			float nDotWo = dot(tri.normal, -wi);

			float t = dot(tri.v1 - attrib.intersection, tri.normal) / dot(wi, tri.normal);
			float3 P = attrib.intersection + t * wi; // intersection in the object space

			float3 tmp0 = tri.v3 - tri.v1;
			float3 tmp1 = tri.v2 - tri.v1;
			float3 tmp2 = P - tri.v1;
			float tmp0dot0 = dot(tmp0, tmp0);
			float tmp0dot1 = dot(tmp0, tmp1);
			float tmp0dot2 = dot(tmp0, tmp2);
			float tmp1dot1 = dot(tmp1, tmp1);
			float tmp1dot2 = dot(tmp1, tmp2);
			float denom = tmp0dot0 * tmp1dot1 - tmp0dot1 * tmp0dot1;

			float u = (tmp1dot1 * tmp0dot2 - tmp0dot1 * tmp1dot2) / denom;
			float v = (tmp0dot0 * tmp1dot2 - tmp0dot1 * tmp0dot2) / denom;

			if (!(0 > u || u > 1 || 0 > v || v > 1 || u + v > 1 || nDotWo == 0.0f || t < 0.001)) {
			    lightDir = wi;
			    light_dist = t;
			    sampled_light_pos = P;
			    shadow_payload.isVisible = true;
			    break;
			}
			else {
			    shadow_payload.isVisible = false;
			}
		    }
		}
		if (!shadow_payload.isVisible) continue;

                Ray shadow_ray = make_Ray(shadow_ray_origin, lightDir, 1, cf.epsilon, light_dist - cf.epsilon);

                rtTrace(root, shadow_ray, shadow_payload);
                float3 n = attrib.normal;

                if (shadow_payload.isVisible) {
                    // rendering equation here: 
		    float3 brdf;
		    float pdf;
		    if (mv.brdf == 0) {
			brdf = (mv.diffuse / M_PIf) +
			    (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                            powf(fmaxf(dot(normalize(reflect(-attrib.wo, attrib.normal)), lightDir), .0f), mv.shininess));
			pdf = ((1.0f - t) * (fmaxf(dot(n, lightDir), .0f) / M_PIf)) + 
				t * ((mv.shininess + 1.0f) / (2.0f * M_PIf)) * 
				powf(fmaxf(dot(r, lightDir), .0f), mv.shininess);
		    }
		    if (mv.brdf == 1) {
			float lightDir_dot_n = dot(lightDir, n);
			float wo_dot_n = dot(attrib.wo, n);
			if (lightDir_dot_n > .0f && wo_dot_n > .0f) {
			    float t = fmaxf(0.25f, s_mean/(s_mean+d_mean));
			    float alpha = mv.roughness;
			    float3 h = normalize(lightDir + attrib.wo); // half angle: 
			    float theta_h = acosf(clamp(dot(h, n),0.0f,1.0f)); // not sure if need to clamp 0
			    // microfacet distribution function, D: 
			    float D = (alpha * alpha) / (M_PIf * powf(cosf(theta_h), 4.0f) *
			           powf((alpha * alpha) + powf(tanf(theta_h), 2.0f), 2.0f));
			    if (powf(cosf(theta_h), 4.0f)*(tanf(theta_h) == 0)) D = 0.0f;
			    // shadow-masking function, G:  
			    float G_1_lightDir = (lightDir_dot_n > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) *
				    powf(tanf(acosf(clamp(dot(lightDir, n),0.0f,1.0f))), 2.0f))) : .0f;
			    float G_1_wo = (wo_dot_n > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) *
				    powf(tanf(acosf(clamp(dot(attrib.wo, n),0.0f,1.0f))), 2.0f))) : .0f;
			    float G = G_1_lightDir * G_1_wo;

			    // fresnel function, F:
			    float3 F = mv.specular + (make_float3(1.0f) - mv.specular) *
				       powf(1.0f - dot(lightDir, h), 5.0f);
			    float3 f_brdf_GGX = (F * G * D) / (4.0f * lightDir_dot_n * wo_dot_n);
			    brdf = (mv.diffuse / M_PIf) + f_brdf_GGX;
			    pdf = fmaxf(((1 - t) * (lightDir_dot_n / M_PIf)) +
			     ((t * D * dot(n, h)) / (4.0f * fmaxf(dot(h, lightDir),0.0f))),0.0f);

		    }
		    else brdf = make_float3(0.0f); // assume f zero otherwise
		}

                    float3 x_prime = sampled_light_pos;
                    float3 x = shadow_ray_origin;
                    float3 n_light = normalize(cross(ab, ac));

                    float R = length(x - x_prime);
		    float pdf_nee;

                    // note: normal should point AWAY from the hitpoint, i.e. dot(n_light, x - x_prime) < 0
                    float G = (1.0f / powf(R, 2.0f)) *
		    fmaxf(dot(n, normalize(x_prime - x)), .0f) *
		    (fmaxf(dot(n_light, normalize(x_prime - x)), .0f));

		    if (cf.NEE == 2 && k < qlights.size() ) {
			pdf_nee = (R*R)/(area*dot(n_light, lightDir))/qlights.size();
			float weight = powf(pdf_nee,2.0f)/(powf(pdf_nee,2.0f) + powf(pdf,2.0f));
			sampled_result = weight * brdf * G * (1.0f/pdf_nee);
			weight = powf(pdf,2.0f)/(powf(pdf_nee,2.0f) + powf(pdf,2.0f));
			sampled_result += weight * brdf * G * (1.0f/pdf);
			sampled_result = sampled_result*100;
			//sampled_result = sampled_result/2;
		    }
		    else if (cf.NEE == 2) {
			pdf_nee = (R*R)/(area*dot(n_light, lightDir))/qlights.size();
			float weight = powf(pdf_nee,2.0f)/(powf(pdf_nee,2.0f) + powf(pdf,2.0f));
			sampled_result += weight * brdf * fmaxf(dot(n,wi),0) * (1.0f/pdf_nee);
			weight = powf(pdf,2.0f)/(powf(pdf_nee,2.0f) + powf(pdf,2.0f));
			sampled_result += weight * brdf * fmaxf(dot(n,wi),0) * (1.0f/pdf);
			//sampled_result = sampled_result/2;
		    }
		    else {
			sampled_result += brdf * G;
		    }

                }
            }
        }
        L_d += qlights[whichLight].color * sampled_result * (area / (float)light_samples);
    }


    
    if (cf.NEE && (payload.depth == 0)) {
        result += L_e;
        payload.radiance = (cf.NEE == 1 ? result + L_d : L_d) * payload.throughput;
    }
    else {
	if (cf.NEE == 2) {
	    result += L_e;
	}
        if (cf.NEE) {
	    result += L_d;
        }
        else {
            result += L_e;
        }
        payload.radiance = result * payload.throughput;
    }

    float q;
    if (cf.RR) {

        q = 1.0f - fmin(fmax(fmax(payload.throughput.x, payload.throughput.y), payload.throughput.z), 1.0f);
        // pick a num from 0 to 1, if less than q, terminate ray
        // i.e. make throughput 0
        if (rnd(payload.seed) < q) {
	    payload.done = true;
	    return;
        }
        else {
            bruh *= (1.0f / (1.0f - q));
        }
    }

    payload.throughput *= bruh;
    payload.origin = attrib.intersection;
    payload.dir = wi;
    payload.depth++;
}
