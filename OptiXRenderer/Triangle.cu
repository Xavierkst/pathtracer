#include "optix.h"
#include "optix_device.h"
#include "Geometries.h"
#include "Payloads.h"

using namespace optix;

rtBuffer<Triangle> triangles; // a buffer of all triangles 

rtDeclareVariable(Ray, ray, rtCurrentRay, );

// Attributes to be passed to material programs 
rtDeclareVariable(Attributes, attrib, attribute attrib, );

// Transfer values into payload variable for color calc in closestHit()
rtDeclareVariable(Payload, payload, rtPayload, );

RT_PROGRAM void intersect(int primIndex)
{
    // Find the intersection of the current ray and triangle
    Triangle tri = triangles[primIndex];
    float t;

    // TODO: implement triangle intersection test here
    
    // get plane normal
    float3 edge1 = tri.vertices[1] - tri.vertices[0]; 
    float3 edge2 = tri.vertices[2] - tri.vertices[0]; 
    // we have normal of triangle, N
    float3 N = normalize(cross(edge1, edge2)); 

    // triangles don't need to worry about transforming
    // ray by M-1 since it is still a triangle (but for spheres, 
    // must consider)

    // find parametric dist t: 
    t = (dot(tri.vertices[0], N) - dot(ray.origin, N)) / dot(ray.direction, N);
    if (t < 0) {
        // triangle is not positive distance to ray, i.e. behind
        return; 
    }

    // check if ray is within the triangle, use barycentric
    // we use cross prod of an edge in the tri and the hit point
    // we know hit point is outside of triangle if dot product of 
    // Normal and orthogonal vect is < 0
    // Note: tri is ACW 
    float3 hitPt = ray.origin + t * ray.direction; 

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
    // cross prod of any 2 edges from above
    float3 hitPtNormal = N;


    // Report intersection (material programs will handle the rest)
    if (rtPotentialIntersection(t))
    {
        // TODO: assign attribute variables here
        // Pass attributes: i.e. materials of the object
        attrib = tri.attributes;

        // Pass hitPt and normal at hitPt into payload 
        // to calculate payload.radiance in closestHit()
        payload.hitPoint = hitPt;
        payload.hitPointNormal = hitPtNormal;

        rtReportIntersection(0);
    }
}

RT_PROGRAM void bound(int primIndex, float result[6])
{
    Triangle tri = triangles[primIndex];

    // TODO: implement triangle bouding box
    result[0] = -1000.f;
    result[1] = -1000.f;
    result[2] = -1000.f;
    result[3] = 1000.f;
    result[4] = 1000.f;
    result[5] = 1000.f;
}