#include "optix.h"
#include "optix_device.h"
#include "Geometries.h"
#include "Payloads.h"

using namespace optix;

rtBuffer<Triangle> triangles; // a buffer of all triangles 

rtDeclareVariable(Ray, ray, rtCurrentRay, );

// Attributes to be passed to material programs 
rtDeclareVariable(Attributes, attrib, attribute attrib, );
rtDeclareVariable(intersectionData, intersectData, attribute intersectData, );

// Transfer values into payload variable for color calc in closestHit()
rtDeclareVariable(Payload, payload, rtPayload, );
//rtDeclareVariable(ShadowPayload, shadowPayload, rtPayload, );

RT_PROGRAM void intersect(int primIndex)
{
    // Find the intersection of the current ray and triangle
    Triangle tri = triangles[primIndex];
    float t;
    float epsilon = 0.001f;

    // TODO: implement triangle intersection test here
    
    // get plane normal
    float3 edge1 = tri.vertices[1] - tri.vertices[0]; 
    float3 edge2 = tri.vertices[2] - tri.vertices[0]; 
    // we have normal of triangle, N
    float3 N = normalize(cross(edge1, edge2)); 

    // triangles don't need to worry about transforming
    // ray by M-1 since it is still a triangle (but for spheres, 
    // must consider)
    float4 temp_ray = tri.transform.inverse() * make_float4(ray.origin, 1);
    float3 ray_orig = make_float3(temp_ray / (float)temp_ray.w);
    //ray_orig = make_float3(ray)
    float3 ray_dir = normalize(make_float3(tri.transform.inverse() * make_float4(ray.direction, 0)));

    //float3 ray_orig = ray.origin;
    //float3 ray_dir = normalize(ray.direction);
    // find parametric dist t: 
    //t = (dot(tri.vertices[0], N) - dot(ray.origin, N)) / dot(ray.direction, N);
    t = (dot(tri.vertices[0], N) - dot(ray_orig, N)) / dot(ray_dir, N);
    if (t < 0) {
        // triangle is not positive distance to ray, i.e. behind
        return; 
    }

    // check if ray is within the triangle, use barycentric
    // we use cross prod of an edge in the tri and the hit point
    // we know hit point is outside of triangle if dot product of 
    // Normal and orthogonal vect is < 0
    // Note: tri is ACW 
    float3 hitPt = ray_orig + t * ray_dir /*+ N*epsilon*/; // account for shadow acne: 
    
    float3 orthogEdge;
    // check 1 (total 3 edges to check) 
    float3 edgeV1V0 = tri.vertices[1] - tri.vertices[0]; 
    float3 edgePV0 = hitPt - tri.vertices[0]; 
    orthogEdge = cross(edgeV1V0, edgePV0); 
    // calc u, v or w here: 

    // check if hitPt inside or outside tri 
    if (dot(N, orthogEdge) < 0)
        return; // hitPt outside

    // Check next edge
    float3 edgeV2V1 = tri.vertices[2] - tri.vertices[1];
    float3 edgePV1 = hitPt - tri.vertices[1];
    orthogEdge = cross(edgeV2V1, edgePV1);

    if (dot(N, orthogEdge) < 0)
        return; 

    // last check: edgeV0V2 and edgePV2
    float3 edgeV0V2 = tri.vertices[0] - tri.vertices[2]; 
    float3 edgePV2 = hitPt - tri.vertices[2];
    orthogEdge = cross(edgeV0V2, edgePV2);

    if (dot(N, orthogEdge) < 0)
        return;

    // made it here, means the hit point is inside the tri!

    // find normal at hitPoint, pass into payload
    // Note: the triangle already transformed during parsing
    //float3 hitPt = ray.origin + t * ray.direction; 
     //cross prod of any 2 edges from above


    // transform normal to worldspace
    float3 hitPtNormal = normalize(make_float3((
        tri.transform.inverse()).transpose() * make_float4(N, 0)));
    // transform hit point to worldspace
    float4 temp_hit = tri.transform * make_float4(hitPt, 1);
    hitPt = make_float3(temp_hit / (float) temp_hit.w) /*+ epsilon * hitPtNormal*/;

    //// obtain parametric distance to hitPoint in worldspace
    t = length(hitPt - ray.origin);

    // compute reflection ray direction
    float3 reflectionDir = normalize(ray.direction - (2.0f * dot(ray.direction, hitPtNormal) * hitPtNormal));
    //float3 reflectionDir = normalize(ray.direction - 2.0f * dot(ray.direction, N) * N);

    // Report intersection (material programs will handle the rest)
    if (rtPotentialIntersection(t))
    {
        // TODO: assign attribute variables here
        // Pass attributes: i.e. materials of the object
        attrib = tri.attributes;
        intersectData.hitPoint = hitPt;
        intersectData.hitPointNormal = /*hitPtNormal*/N;
        intersectData.reflectDir = reflectionDir;
        intersectData.rayDir = ray.direction;
        intersectData.rayOrig = ray.origin;
        // Pass hitPt and normal at hitPt into payload 
        // to calculate payload.radiance in closestHit()
        //payload.hitPoint = hitPt;
        //payload.hitPointNormal = hitPtNormal;
        //payload.dir = reflectionDir;

        //shadowPayload.isVisible = true;

        rtReportIntersection(0);
    }
}

RT_PROGRAM void bound(int primIndex, float result[6])
{
    Triangle tri = triangles[primIndex];

    result[0] = -1000.f;
    result[1] = -1000.f;
    result[2] = -1000.f;
    result[3] = 1000.f;
    result[4] = 1000.f;
    result[5] = 1000.f;

    // TODO: implement triangle bouding box
    if (tri.vertices[0].x > result[0]) {
        result[0] = tri.vertices[0].x;
    }
    
    if (tri.vertices[1].x > result[0]) {
        result[0] = tri.vertices[1].x;
    }

    if (tri.vertices[2].x > result[0]) {
        result[0] = tri.vertices[2].x;
    }

    // y 
    if (tri.vertices[0].y > result[0]) {
        result[0] = tri.vertices[0].y;
    }
    
    if (tri.vertices[1].y > result[0]) {
        result[0] = tri.vertices[1].y;
    }

    if (tri.vertices[2].y > result[0]) {
        result[0] = tri.vertices[2].y;
    }
    
    // z
    if (tri.vertices[0].z > result[0]) {
        result[0] = tri.vertices[0].z;
    }
    
    if (tri.vertices[1].z > result[0]) {
        result[0] = tri.vertices[1].z;
    }

    if (tri.vertices[2].z > result[0]) {
        result[0] = tri.vertices[2].z;
    }

}