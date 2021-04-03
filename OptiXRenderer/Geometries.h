#pragma once

#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>

 /* 
 * The attributes defining what material the geoms are 
 * made of
 */
struct Attributes
{
    // TODO: define the attributes structure
    float shininess; 
    optix::float3 diffuse;
    optix::float3 specular;
    optix::float3 emission;
};

/**
 * Structures describing different geometries should be defined here.
 */
struct Triangle
{
    // TODO: define the triangle structure
    //int maxverts; 
    //int maxvertnorms;
    optix::float3 vertices[3];
    optix::float3 normals[3];
    //int triIndices[3];
    //int normIndices[3];

    // include triangle materials here: 
    optix::float3 ambient;
    Attributes attributes;

    // contains the transf -- used to inverse transform
    // the ray later (just peek the top of stack, and assign to transform)
    optix::Matrix4x4 transform;

    Triangle() {
        for (int i = 0; i < 3; i++) {
            vertices[i] = optix::make_float3(.0f);
            normals[i] = optix::make_float3(.0f);
            //triIndices[i] = 0;
            //normIndices[i] = 0;
        }

        ambient = optix::make_float3(.0f);
        attributes.shininess = .0f;
        attributes.diffuse = optix::make_float3(.0f);
        attributes.specular = optix::make_float3(.0f);
        attributes.emission = optix::make_float3(.0f);
    }

    Triangle(optix::float3 ) {
        for (int i = 0; i < 3; i++) {
            vertices[i] = optix::make_float3(.0f);
            normals[i] = optix::make_float3(.0f);
            //triIndices[i] = 0;
            //normIndices[i] = 0;
        }

        ambient = optix::make_float3(.0f);
        attributes.shininess = .0f;
        attributes.diffuse = optix::make_float3(.0f);
        attributes.specular = optix::make_float3(.0f);
        attributes.emission = optix::make_float3(.0f);
    }
};

struct Sphere
{
    // TODO: define the sphere structure
    optix::float3 position;
    float radius;

    // include Sphere materials here:
    optix::float3 ambient;
    Attributes attributes;

    // contains the transf -- used to inverse transform
    // the ray later
    optix::Matrix4x4 transform;

    Sphere() {
        position = optix::make_float3(.0f);
        radius = .0f;

        ambient = optix::make_float3(.0f);
        attributes.shininess = .0f;
        attributes.diffuse = optix::make_float3(.0f);
        attributes.specular = optix::make_float3(.0f);
        attributes.emission = optix::make_float3(.0f);
    }
};

