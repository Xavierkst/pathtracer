#pragma once

#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>

/**
 * Structures describing different light sources should be defined here.
 */

struct Attenuation {
    float constant;
    float linear; 
    float quadratic;

    Attenuation() {
        constant = 1.0f;
        linear = .0f;
        quadratic = .0f;
    }
};

/*
 * All lights can take in an attenutation value
 * else, no attenuation is 1, 0, 0 (const, lin, quad) by default
 */
struct PointLight
{
    // TODO: define the point light structure
    // needs a position, light color
    optix::float3 light_pos;
    optix::float3 light_color;
    
    Attenuation attenuation; 

    PointLight(optix::float3 position, optix::float3 color,
        Attenuation atten) {
        light_pos = position;
        light_color = color;
        attenuation.constant = atten.constant;
        attenuation.linear = atten.linear;
        attenuation.quadratic = atten.quadratic;
    }
};

struct DirectionalLight
{
    // TODO: define the directional light structure
    // needs a light-direction, light color
    optix::float3 light_dir;
    optix::float3 light_color;

    DirectionalLight(optix::float3 dir, optix::float3 color) {
        light_dir = optix::normalize(dir);
        light_color = color;
    }
};