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
rtDeclareVariable(uint, next_event_est, , );
rtDeclareVariable(uint, sampling_method, , );

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
        float3 sampled_result = make_float3(.0f);
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

RT_PROGRAM void pathTracer() {

    MaterialValue mv = attrib.mv;
    Config cf = config[0];

    float3 n = normalize(attrib.normal);
    float3 L_e = mv.emission;
    float3 result = make_float3(.0f);
    float3 L_d = make_float3(.0f);
    float3 K_s = mv.specular;
    float3 K_d = mv.diffuse;
    float3 r = normalize(reflect(-attrib.wo, attrib.normal));

    float exp_beta = 2.0f;

    // variables for direct lighting
    int root_light_samples = (int)sqrtf(light_samples);
    float pdf_NEE = .0f;

    // When next event estimation is ON:
    // if an indir ray ever strikes light source (and it is NOT the first ray cast)
    // ray should be terminated
    if (cf.next_event_est && attrib.objType == LIGHT && payload.depth != 0) {
        payload.depth = cf.maxDepth;
        payload.done = true;
        payload.radiance = result;
        return;
    }

    // Add indirect lighting here:
    // generate randomize ray direction w_i
    float zeta_0 = rnd(payload.seed);
    float zeta_1 = rnd(payload.seed);
    float zeta_2 = rnd(payload.seed);
    float K_s_avg = (K_s.x + K_s.y + K_s.z) / 3.0f;
    float K_d_avg = (K_d.x + K_d.y + K_d.z) / 3.0f;
    //rtPrintf("ks %f\n", K_s_avg);
    //rtPrintf("kd %f\n", K_d_avg);
    float t = .0f;
    // make sure denom not 0, else check if phong or ggx
    // and set t accordingly
    if ((K_s_avg == .0f) && (K_d_avg == .0f)) {
        //rtPrintf("heeeere\n");
        // if denom is 0, t is 0 when brdf is modified phong, and 1 for ggx
        if (mv.brdf_type == MOD_PHONG)
            t = .0f;
        else if (mv.brdf_type == GGX)
            t = 1.0f;
    }
    else {
        // denom is not zero
        //rtPrintf("OR heeeere\n");
        //rtPrintf("ks %f\n", K_s_avg);
        //rtPrintf("kd %f\n", K_d_avg);
        float k_val = K_s_avg / (K_s_avg + K_d_avg);
        if (mv.brdf_type == MOD_PHONG) {
            t = k_val;
            //rtPrintf("k value is: %f\n", k_val);
            //if (t != .0f) 
            //    rtPrintf("t value is: %f\n", t);

        }
        else {
            t = fmaxf(.25f, k_val);
            //rtPrintf("t value is: %f\n", t);
        }
    }

    float theta = .0f;
    float phi = .0f;
    switch (sampling_method) {
    case HEMISPHERE_SAMPLING:
        //rtPrintf("hemisphere here\n");
        phi = 2.0f * M_PIf * zeta_2;
        theta = acosf(zeta_1);
        break;
    case COSINE_SAMPLING:
        //rtPrintf("cosine here\n");
        phi = 2.0f * M_PIf * zeta_2;
        theta = acosf(sqrtf(zeta_1));
        break;
    case BRDF_SAMPLING:
        //rtPrintf("brdf here\n");
        // phi remains the same for either specular or diffuse pdf
        phi = 2.0f * M_PIf * zeta_2;
        if (zeta_0 > t)
            theta = acosf(sqrtf(zeta_1)); // theta_diffuse
        else
            theta = acosf(powf(zeta_1, (1.0f / (mv.shininess + 1.0f)))); // theta_specular
        break;
    default:
    }

    // qn: why rotate s wrt the z-axis? and not the y-axis?
    float3 sample_s = make_float3(cosf(phi) * sinf(theta), sinf(phi) * sinf(theta), cosf(theta));

    // generate coordinate frame at the intersect point
    //float3 n = normalize(attrib.normal);
    float3 w = n;
    // if specular brdf chosen, center sample_s at reflection vector r
    if ((zeta_0 <= t) && (sampling_method == BRDF_SAMPLING)) w = r;

    //else w = normalize(n);

    float3 a = make_float3(.0f, 1.0f, .0f);
    // incase a and w are closely aligned, swap a out for 
    // a diff arbitrary vector <1,0,0> instead of <0,1,0>
    if (1.0f - fabsf(dot(a, w)) <= 1.0f) {
        a = make_float3(1.0f, .0f, .0f);
    }

    float3 u = normalize(cross(a, w));
    float3 v = normalize(cross(w, u)); // i dont think need to normalize

    float theta_h_sample = atanf((mv.roughness * sqrtf(zeta_1)) / sqrtf(1.0f - zeta_1));
    float phi_h_sample = 2.0f * M_PIf * zeta_2;
    float3 h_sample = make_float3(cosf(phi_h_sample) * sinf(theta_h_sample),
        sinf(phi_h_sample) * sinf(theta_h_sample),
        cosf(theta_h_sample));

    float3 w_h = n;
    float3 a_h = make_float3(.0f, 1.0f, .0f);
    if (1.0f - fabsf(dot(a_h, w_h)) <= 1.0f) {
        a_h = make_float3(1.0f, .0f, .0f);
    }

    float3 u_h = normalize(cross(a_h, w_h));
    float3 v_h = normalize(cross(w_h, u_h));

    // rotate h:
    float3 w_i = make_float3(.0f);
    // get randomized new ray dir -- choose a sampling method
    if (mv.brdf_type == MOD_PHONG)
        w_i = normalize((sample_s.x * u + sample_s.y * v + sample_s.z * w));
    else {
        h_sample = (h_sample.x * u_h + h_sample.y * v_h + h_sample.z * w_h);
        w_i = normalize(reflect(-attrib.wo, h_sample));
    }

    if (cf.next_event_est) {
        // Add direct lighting here:
        for (int k = 0; k < qlights.size(); ++k) {
            float3 sampled_result = make_float3(.0f);
            float pdf_lights_k = .0f;
            float pdf_brdf = .0f;
            float3 w_i_dir = make_float3(.0f);
            float3 h = make_float3(.0f);
            float D = .0f;
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
                        //float D = .0f;

                        // rendering equation here: 
                        //float3 w_i = sampled_light_pos;
                        w_i_dir = normalize(sampled_light_pos - shadow_ray_origin);

                        float3 f_brdf = make_float3(.0f);
                        if (mv.brdf_type == MOD_PHONG) {
                            f_brdf = (mv.diffuse / M_PIf) +
                                (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                                    powf(fmaxf(dot(r, w_i_dir), .0f), mv.shininess));
                        }
                        else {
                            float wi_dot_n_dir = dot(w_i_dir, n);
                            float wo_dot_n_dir = dot(attrib.wo, n);
                            if (wi_dot_n_dir > .0f && wo_dot_n_dir > .0f) {
                                float alpha = mv.roughness;
                                float3 h = normalize(w_i_dir + attrib.wo); // half angle: 
                                float theta_h = acosf(dot(h, n)); // not sure if need to clamp 0
                                float cos_theta_h_4 = powf(cosf(theta_h), 4.0f);
                                float alpha_tan_theta_h_sq = (alpha * alpha) + powf(tanf(theta_h), 2.0f);
                                // make sure denom not 0, else set D to 0
                                if (cos_theta_h_4 * alpha_tan_theta_h_sq != .0f) D = (alpha * alpha) / (M_PIf * cos_theta_h_4 * powf(alpha_tan_theta_h_sq, 2.0f));

                                // microfacet distribution function, D: 
                                //float D = (alpha * alpha) / (M_PIf * powf(cosf(theta_h), 4.0f) * powf((alpha * alpha) + powf(tanf(theta_h), 2.0f), 2.0f));

                                // shadow-masking function, G:  
                                float G_1_wi = (wi_dot_n_dir > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) * powf(tanf(acosf(dot(w_i_dir, n))), 2.0f))) : .0f;
                                float G_1_wo = (wo_dot_n_dir > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) * powf(tanf(acosf(dot(attrib.wo, n))), 2.0f))) : .0f;
                                float G = G_1_wi * G_1_wo;
                                // fresnel function, F:
                                float3 F = K_s + (1.0f - K_s) * powf(1.0f - dot(w_i_dir, h), 5.0f);
                                float3 f_brdf_GGX = (F * G * D) / (4.0f * wi_dot_n_dir * wo_dot_n_dir);

                                f_brdf = (K_d / M_PIf) + f_brdf_GGX;
                            }
                        }

                        float3 x_prime = sampled_light_pos;
                        float3 x = shadow_ray_origin;
                        float3 n = attrib.normal;
                        float3 n_light = normalize(cross(ab, ac));

                        float R = length(x - x_prime);

                        // note: normal should point AWAY from the hitpoint, i.e. dot(n_light, x - x_prime) < 0
                        float G = (1.0f / powf(R, 2.0f)) * fmaxf(dot(n, normalize(x_prime - x)), .0f) *
                            (fmaxf(dot(n_light, normalize(x_prime - x)), .0f));

                        sampled_result += f_brdf * G;
                        //pdf_lights_k += (powf(R, 2.0f) / (area * fabsf(dot(n, w_i_dir)))); 
                    }
                }
            }
            // calculate weight Wi for this given w_i generated
            pdf_NEE = pdf_lights_k /** (1.0f / (float) qlights.size())*/;

            // calculate pdf_brdf 
            pdf_brdf = ((1.0f - t) * fmaxf(dot(n, w_i_dir), .0f) / M_PIf) + ((t * D * dot(n, h)) / (4.0f * dot(h, w_i_dir)));
         
            // calc power heuristic: 
            float pdf_denom_sum = powf(pdf_NEE, exp_beta) + powf(pdf_brdf, exp_beta);
            float pdf_numerator = powf(pdf_NEE, exp_beta);
            float weight_i = pdf_numerator / pdf_denom_sum;
        
            if (cf.next_event_est == MIS) {
                //rtPrintf("here");
                L_d += qlights[k].color * sampled_result * (area / (float)light_samples) * (1.0f / pdf_NEE) * weight_i;
            }
            else {
                // divide brdf by the pdf here
                L_d += qlights[k].color * sampled_result * (area / (float)light_samples);
            }
        }
    }

    //// Add indirect lighting here:
    //// generate randomize ray direction w_i
    //float zeta_0 = rnd(payload.seed);
    //float zeta_1 = rnd(payload.seed);
    //float zeta_2 = rnd(payload.seed);
    //float K_s_avg = (K_s.x + K_s.y + K_s.z) / 3.0f;
    //float K_d_avg = (K_d.x + K_d.y + K_d.z) / 3.0f;
    ////rtPrintf("ks %f\n", K_s_avg);
    ////rtPrintf("kd %f\n", K_d_avg);
    //float t = .0f;
    //// make sure denom not 0, else check if phong or ggx
    //// and set t accordingly
    //if ((K_s_avg == .0f) && (K_d_avg == .0f)) {
    //    //rtPrintf("heeeere\n");
    //    // if denom is 0, t is 0 when brdf is modified phong, and 1 for ggx
    //    if (mv.brdf_type == MOD_PHONG) 
    //        t = .0f;
    //    else if (mv.brdf_type == GGX) 
    //        t = 1.0f;
    //}
    //else {
    //    // denom is not zero
    //    //rtPrintf("OR heeeere\n");
    //    //rtPrintf("ks %f\n", K_s_avg);
    //    //rtPrintf("kd %f\n", K_d_avg);
    //    float k_val = K_s_avg / (K_s_avg + K_d_avg);
    //    if (mv.brdf_type == MOD_PHONG) {
    //        t = k_val;
    //        //rtPrintf("k value is: %f\n", k_val);
    //        //if (t != .0f) 
    //        //    rtPrintf("t value is: %f\n", t);

    //    }
    //    else {
    //        t = fmaxf(.25f, k_val);
    //        //rtPrintf("t value is: %f\n", t);
    //    }
    //}

    //float theta = .0f;
    //float phi = .0f;
    //switch (sampling_method) {
    //    case HEMISPHERE_SAMPLING: 
    //        //rtPrintf("hemisphere here\n");
    //        phi = 2.0f * M_PIf * zeta_2;
    //        theta = acosf(zeta_1);
    //        break;
    //    case COSINE_SAMPLING: 
    //        //rtPrintf("cosine here\n");
    //        phi = 2.0f * M_PIf * zeta_2;
    //        theta = acosf(sqrtf(zeta_1));
    //        break;
    //    case BRDF_SAMPLING: 
    //        //rtPrintf("brdf here\n");
    //        // phi remains the same for either specular or diffuse pdf
    //        phi = 2.0f * M_PIf * zeta_2;
    //        if (zeta_0 > t) 
    //            theta = acosf(sqrtf(zeta_1)); // theta_diffuse
    //        else 
    //            theta = acosf(powf(zeta_1, (1.0f / (mv.shininess + 1.0f)))); // theta_specular
    //        break;
    //    default: 
    //        break;
    //}

    //// qn: why rotate s wrt the z-axis? and not the y-axis?
    //float3 sample_s = make_float3(cosf(phi) * sinf(theta), sinf(phi) * sinf(theta), cosf(theta));
    //
    //// generate coordinate frame at the intersect point
    //float3 w = n;
    //// if specular brdf chosen, center sample_s at reflection vector r
    //if ((zeta_0 <= t) && (sampling_method == BRDF_SAMPLING)) w = r;

    ////else w = normalize(n);

    //float3 a = make_float3(.0f, 1.0f, .0f);
    //// incase a and w are closely aligned, swap a out for 
    //// a diff arbitrary vector <1,0,0> instead of <0,1,0>
    //if (1.0f - fabsf(dot(a, w)) <= 1.0f) {
    //    a = make_float3(1.0f, .0f, .0f);
    //}

    //float3 u = normalize(cross(a, w));
    //float3 v = normalize(cross(w, u)); // i dont think need to normalize

    //float theta_h_sample = atanf((mv.roughness * sqrtf(zeta_1)) / sqrtf(1.0f - zeta_1));
    //float phi_h_sample = 2.0f * M_PIf * zeta_2;
    //float3 h_sample = make_float3(cosf(phi_h_sample) * sinf(theta_h_sample), 
    //    sinf(phi_h_sample)*sinf(theta_h_sample), 
    //        cosf(theta_h_sample));

    //float3 w_h = n; 
    //float3 a_h = make_float3(.0f, 1.0f, .0f); 
    //if (1.0f - fabsf(dot(a_h, w_h)) <= 1.0f) {
    //    a_h = make_float3(1.0f, .0f, .0f);
    //}
    //
    //float3 u_h = normalize(cross(a_h, w_h)); 
    //float3 v_h = normalize(cross(w_h, u_h));

    //// rotate h:
    //float3 w_i = make_float3(.0f); 
    //// get randomized new ray dir -- choose a sampling method
    //if (mv.brdf_type == MOD_PHONG) 
    //    w_i = normalize((sample_s.x * u + sample_s.y * v + sample_s.z * w));
    //else {
    //    h_sample = (h_sample.x * u_h + h_sample.y * v_h + h_sample.z * w_h);
    //    w_i = normalize(reflect(-attrib.wo, h_sample));
    //}



    // calculate the summation of all pdfs here (i.e. pdf_nee + pdf_brdf)
    // using h, w_i, 
    
    // the BRDF 
    float3 f_brdf = make_float3(0.0f);
    float pdf = 1.0f;
    float3 addon_throughput = make_float3(.0f);

    switch (sampling_method) {
    case HEMISPHERE_SAMPLING:
        //rtPrintf("hemisphere here");
        f_brdf = (mv.diffuse / M_PIf) +
            (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                powf(fmaxf(dot(r, w_i), .0f), mv.shininess));
        pdf = 1.0f / (2.0f * M_PIf);
        addon_throughput = (f_brdf * fmaxf(dot(n, w_i), .0f) * (1.0f / pdf));
        break;
    case COSINE_SAMPLING:
        //rtPrintf("cosine here");
        f_brdf = (mv.diffuse / M_PIf) +
            (mv.specular * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                powf(fmaxf(dot(r, w_i), .0f), mv.shininess));
        pdf = fmaxf(dot(n, w_i), .0f) / (M_PIf);
        addon_throughput = (f_brdf)*fmaxf(dot(n, w_i), .0f) * (1.0f / pdf);
        break;
    case BRDF_SAMPLING:
        float pdf_NEE_brdf = .0f;
        //rtPrintf("brdf here");
        // check the material whether to use mod-phong or GGX brdf
        if (mv.brdf_type == MOD_PHONG) {
            f_brdf = (K_d / M_PIf) +
                (K_s * ((mv.shininess + 2.0f) / (2.0f * M_PIf)) *
                    powf(fmaxf(dot(r, w_i), .0f), mv.shininess));

            pdf = ((1.0f - t) * (dot(n, w_i) / M_PIf)) +
                t * ((mv.shininess + 1.0f) / (2.0f * M_PIf)) *
                powf(fmaxf(dot(r, w_i), .0f), mv.shininess);
        }
        else {
            // construct GGX BRDF: 
            float pdf_lights_k = .0f;

            float wi_dot_n = dot(w_i, n);
            float wo_dot_n = dot(attrib.wo, n);
            if (wi_dot_n > .0f && wo_dot_n > .0f) {
                float alpha = mv.roughness;
                float3 h = normalize(w_i + attrib.wo); // half angle: 
                float theta_h = acosf(dot(h, n)); // not sure if need to clamp 0
                // microfacet distribution function, D: 
                float cos_theta_h_4 = powf(cosf(theta_h), 4.0f);
                float alpha_tan_theta_h_sq = (alpha * alpha) + powf(tanf(theta_h), 2.0f);
                float D = .0f;
                // make sure denom not 0, else set D to 0
                if (cos_theta_h_4 * alpha_tan_theta_h_sq != .0f) D = (alpha * alpha) / (M_PIf * cos_theta_h_4 * powf(alpha_tan_theta_h_sq, 2.0f));
                //D = (alpha * alpha) / (M_PIf * powf(cosf(theta_h), 4.0f) * powf(alpha_tan_theta_h_sq, 2.0f));
            // shadow-masking function, G:  
                float G_1_wi = (wi_dot_n > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) * powf(tanf(acosf(dot(w_i, n))), 2.0f))) : .0f;
                float G_1_wo = (wo_dot_n > .0f) ? 2.0f / (1.0f + sqrtf(1.0f + (alpha * alpha) * powf(tan(acosf(dot(attrib.wo, n))), 2.0f))) : .0f;
                float G = G_1_wi * G_1_wo;
                // fresnel function, F:
                float3 F = K_s + (1.0f - K_s) * powf(1.0f - fmaxf(dot(w_i, h), .0f), 5.0f);
                float3 f_brdf_GGX = (F * G * D) / (4.0f * wi_dot_n * wo_dot_n);
                f_brdf = (K_d / M_PIf) + f_brdf_GGX;
                pdf = ((1.0f - t) * fmaxf(dot(n, w_i), .0f) / M_PIf) + ((t * D * dot(n, h)) / (4.0f * dot(h, w_i)));

                for (int k = 0; k < qlights.size(); ++k) {
                    //float3 sampled_result = make_float3(.0f);
                    float pdf_brdf = .0f;
                    // Compute direct lighting equation for w_i_k ray, for k = 1 to N*N
                    float3 a = qlights[k].tri1.v1;
                    float3 b = qlights[k].tri1.v2;
                    float3 c = qlights[k].tri2.v3;
                    float3 d = qlights[k].tri2.v2;

                    float3 ac = c - a;
                    float3 ab = b - a;
                    float area = length(cross(ab, ac));


                    float3 shadow_ray_origin = attrib.intersection /*+ attrib.normal * cf.epsilon*/;
                    float3 shadow_ray_dir = w_i;
                    // trace the ray and see if it hits a light source
                    Ray shadow_ray = make_Ray(shadow_ray_origin, shadow_ray_dir, 1, cf.epsilon, RT_DEFAULT_MAX);
                    //float light_dist = length(sampled_light_pos - shadow_ray_origin);

                    ShadowPayload shadow_payload;
                    shadow_payload.isVisible = true;
                    rtTrace(root, shadow_ray, shadow_payload);

                    if (!shadow_payload.isVisible && shadow_payload.objType == LIGHT) {
                        float3 x_prime = shadow_payload.intersectPt;
                        float3 x = shadow_ray_origin;
                        float3 n_light = normalize(cross(ab, ac));
                        float R = length(x - x_prime);

                        pdf_lights_k += (powf(R, 2.0f) / (area * fabsf(dot(n, w_i))));
                    }
                    /*else {
                        pdf_lights_k += .0f;
                    }*/
                }

            }
            else f_brdf = make_float3(.0f); // assume f zero otherwise
            pdf_NEE_brdf = /*(1.0f / qlights.size()) **/ pdf_lights_k;
        }

        float pdf_denom_sum = powf(pdf_NEE_brdf, exp_beta) + powf(pdf, exp_beta);
        float pdf_numerator = powf(pdf, exp_beta);
        float weight_i_brdf = pdf_numerator / pdf_denom_sum;

        if (cf.next_event_est == MIS) {
            addon_throughput = (f_brdf * fmaxf(dot(n, w_i), .0f) * (1.0f / pdf) * weight_i_brdf);
        }
        else {
            addon_throughput = (f_brdf * fmaxf(dot(n, w_i), .0f) * (1.0f / pdf));
        }
        break;
        default:
    }


    // Check if its first intersected surface
    if (cf.next_event_est && (payload.depth == 0)) {
        result += L_e;
        payload.radiance = (result + L_d) * payload.throughput;
    }
    // on the last bounce, we return only emission term
    // NEE true (1) or false (0) 
    // stop recursion at depth D-1 if NEE is true
    //else if (payload.depth == (cf.maxDepth /*- 1 - cf.next_event_est*/)) {
    //    //rtPrintf("is this ever called?");
    //    result += L_e;
    //    payload.radiance = result;
    //    payload.done = true;
    //}
    else {
        if (cf.next_event_est) {
            //if (attrib.objType == LIGHT) result = make_float3(.0f);
            result += L_d;
        }
        else {
            result += L_e;
        }
        payload.radiance = result * payload.throughput;
    }

    // Apply russian roulette to determine whether to terminate ray
    // or boost its contribution
    float q;
    if (cf.russian_roul) {
        q = 1.0f - fmin(fmax(fmax(payload.throughput.x, payload.throughput.y), payload.throughput.z), 1.0f);
        // pick a num from 0 to 1, if less than q, terminate ray
        // i.e. make throughput 0
        if (rnd(payload.seed) < q) {
            addon_throughput *= make_float3(.0f);
        }
        else {
            float thru_put_boost = (1.0f / (1.0f - q));
            addon_throughput *= thru_put_boost;
        }
    }
    payload.throughput *= addon_throughput;
    payload.origin = attrib.intersection;
    payload.dir = w_i;

    payload.depth++;
}

