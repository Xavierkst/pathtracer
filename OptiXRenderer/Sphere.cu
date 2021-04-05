#include <optix.h>
#include <optix_device.h>
#include "Geometries.h"
#include "Payloads.h"


using namespace optix;

rtBuffer<Sphere> spheres; // a buffer of all spheres

rtDeclareVariable(Ray, ray, rtCurrentRay, );

// Attributes to be passed to material programs 
rtDeclareVariable(Attributes, attrib, attribute attrib, );

// pass to payload: hitPt, Normal at hitPt
rtDeclareVariable(Payload, payload, rtPayload, );

RT_PROGRAM void intersect(int primIndex)
{
    // Find the intersection of the current ray and sphere
    Sphere sphere = spheres[primIndex];
    float t, t_sol1, t_sol2;
    float a, b, c;
    // TODO: implement sphere intersection test here
    // there are 4 outcomes for the ray here
    // completely miss. ray pass intersect 2 pts. 
    // ray skims the surface. 
    // ray intersects only 1 pt, i.e. ray starts from inside

    // apply the inverse transform to ray here before calculating 
    // sphere intersection
    float3 ray_orig = make_float3(sphere.transform.inverse() * make_float4(ray.origin, 1));
    float3 ray_dir = normalize(make_float3(sphere.transform.inverse() * make_float4(ray.direction, 0)));

    // use quadratic equation and find solutions for t 
    optix::float3 rad = ray_orig - sphere.position; 
    a = dot(ray_dir, ray_dir); 
    b = 2.0f * dot(ray_dir, rad);
    c = dot(rad, rad) - (sphere.radius * sphere.radius);

    // find discriminant i.e. the value inside sqrt (b^2 - 4ac) 
    float disc = (b * b) - (4 * a * c);

    // if discrim is negative, no intersection
    if (disc < .0f) return; 

    disc = sqrtf(disc); 
    // both values are legal
    t_sol1 = (-1.0f * b + disc ) / (2.0f * a);
    t_sol2 = (-1.0f * b - disc ) / (2.0f * a);
    
    // find the smaller of the 2
    if (t_sol2 < t_sol1) {
        float t_temp = t_sol1;
        t_sol1 = t_sol2;
        t_sol2 = t_temp;
    }
    // now check to make sure we chose a positive param-distance
    if (t_sol1 < .0f) t_sol1 = t_sol2;
    // both solutions are negative
    if (t_sol1 < .0f) return; 

    t = t_sol1;

    // find normal at hitPoint, pass into payload
    float3 hitPt = ray_orig + t * ray_dir; 
    float3 hitPtNormal = normalize(hitPt - sphere.position);

    // Report intersection (material programs will handle the rest)
    if (rtPotentialIntersection(t))
    {
        // TODO: assign attribute variables here
        // Pass attributes
        attrib = sphere.attributes;
        //rtPrintf("made it! %f", hitPt.x);
        // Pass hitPt and normal at hitPt into payload 
        // to calculate payload.radiance in closestHit()

        // transform hitPt back: Mp 
        payload.hitPoint = 
            make_float3(sphere.transform * make_float4(hitPt, 1)); // applying Mp 
        // transform normal at hitPoint back (M-1)^T
        payload.hitPointNormal = normalize(
            make_float3((sphere.transform.inverse()).transpose() * 
                make_float4(hitPtNormal, 0)));

        rtReportIntersection(0);
    }
}

RT_PROGRAM void bound(int primIndex, float result[6])
{
    Sphere sphere = spheres[primIndex];

    // TODO: implement sphere bouding box
    result[0] = -1000.f;
    result[1] = -1000.f;
    result[2] = -1000.f;
    result[3] = 1000.f;
    result[4] = 1000.f;
    result[5] = 1000.f;
}