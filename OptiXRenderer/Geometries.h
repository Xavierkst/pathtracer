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
    optix::float3 ambient;
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

    
    //optix::float3 ambient; // include triangle materials here: 
    Attributes attributes;

    // contains the transf -- used to inverse transform
    // the ray later (just peek the top of stack, and assign to transform)
    optix::Matrix4x4 transform; // used for inv. transf of ray

    // constructor for tris w/ vertex, but no normals
    Triangle(optix::float3 v1, optix::float3 v2, optix::float3 v3,
                optix::float3 amb, float shine, optix::float3 diff, 
                optix::float3 spec, optix::float3 emiss) {
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;

        normals[0] = optix::make_float3(.0f);
        normals[1] = optix::make_float3(.0f);
        normals[2] = optix::make_float3(.0f);

        attributes.ambient = amb;
        attributes.shininess = shine;
        attributes.diffuse = diff;
        attributes.specular = spec;
        attributes.emission = emiss;
    }

    // constructor for tris w/ vertex + normals
    Triangle(optix::float3 v1, optix::float3 v2, optix::float3 v3,
                optix::float3 n1,  optix::float3 n2, optix::float3 n3, 
                optix::float3 amb, float shine, optix::float3 diff, 
                optix::float3 spec, optix::float3 emiss) {
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;

        normals[0] = n1;
        normals[1] = n2;
        normals[2] = n3;

        attributes.ambient = amb;
        attributes.shininess = shine;
        attributes.diffuse = diff;
        attributes.specular = spec;
        attributes.emission = emiss;
    }

    Triangle() {
        for (int i = 0; i < 3; i++) {
            vertices[i] = optix::make_float3(.0f);
            normals[i] = optix::make_float3(.0f);
        }

        attributes.ambient = optix::make_float3(.2f);
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

    Sphere(optix::float3 pos, float r,
        optix::float3 amb, float shine, optix::float3 diff,
        optix::float3 spec, optix::float3 emiss) {
        
        position = pos;
        radius = r;
        ambient = amb;

        attributes.shininess = shine;
        attributes.diffuse = diff;
        attributes.specular = spec;
        attributes.emission = emiss;
    }

    Sphere() {
        position = optix::make_float3(.0f);
        radius = .0f;

        ambient = optix::make_float3(.2f);
        attributes.shininess = .0f;
        attributes.diffuse = optix::make_float3(.0f);
        attributes.specular = optix::make_float3(.0f);
        attributes.emission = optix::make_float3(.0f);
    }
};

