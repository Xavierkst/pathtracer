#include <optix.h>
#include <optix_device.h>
#include "Geometries.h"
#include "Payloads.h"


using namespace optix;

rtBuffer<Sphere> spheres; // a buffer of all spheres

rtDeclareVariable(Ray, ray, rtCurrentRay, );

// Attributes to be passed to material programs 
rtDeclareVariable(Attributes, attrib, attribute attrib, );
rtDeclareVariable(intersectionData, intersectData, attribute intersectData, );

// pass to payload: hitPt, Normal at hitPt
rtDeclareVariable(Payload, payload, rtPayload, );
//rtDeclareVariable(ShadowPayload, shadowPayload, rtPayload, );

RT_PROGRAM void intersect(int primIndex)
{
    // Find the intersection of the current ray and sphere
    Sphere sphere = spheres[primIndex];
    float t, t_sol1, t_sol2;
    float a, b, c;
    float epsilon = 0.001f;

    // TODO: implement sphere intersection test here
    // there are 4 outcomes for the ray here
    // completely miss. ray pass intersect 2 pts. 
    // ray skims the surface. 
    // ray intersects only 1 pt, i.e. ray starts from inside

    // apply the inverse transform to ray here before calculating 
    // sphere intersection
    float4 temp_ray_orig = sphere.transform.inverse() * make_float4(ray.origin, 1);
    float3 ray_orig = make_float3(temp_ray_orig / (float)temp_ray_orig.w);
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

    // obtain normal
    //float3 hitPtNormal = normalize( hitPt - make_float3(sphere.transform.inverse() * make_float4(sphere.position, 1.0f)));
    float3 hitPtNormal = normalize(hitPt - sphere.position);
    //hitPt += epsilon * hitPtNormal;
    // transform normal to worldspace
    hitPtNormal = normalize(make_float3(
        (sphere.transform.inverse()).transpose() * make_float4(hitPtNormal, 0)));

    // transform hitPoint back into worldspace from local space
    float4 temp_pt = sphere.transform * make_float4(hitPt, 1);
    hitPt = make_float3(temp_pt / (float)temp_pt.w)/* + epsilon * hitPtNormal*/;/*+epsilon * hitPtNormal;*/

    // find distance t:  
    t = length(hitPt - ray.origin);
    
    // Account for shadow Acne??  
    //hitPt += hitPtNormal * epsilon;

    // compute reflection ray direction in world space
    //float3 reflectionDir = normalize(ray_dir - 2.0f * dot(ray_dir, hitPtNormal) * hitPtNormal);
    float3 reflectionDir = normalize(ray.direction - (2.0f * dot( ray.direction, hitPtNormal) * hitPtNormal));
    //rtPrintf("%f, %f, %f", reflectionDir.x, reflectionDir.y, reflectionDir.z);

    // Report intersection (material programs will handle the rest)
    if (rtPotentialIntersection(t))
    {
        // TODO: assign attribute variables here
        // Pass attributes
        attrib = sphere.attributes;       
        //rtPrintf("sphere spec: %f %f %f\n", attrib.specular.x, attrib.specular.y, attrib.specular.z);
        intersectData.hitPoint = hitPt;
        
        intersectData.hitPointNormal = hitPtNormal;
        intersectData.reflectDir = reflectionDir;
        intersectData.rayDir = ray.direction;
        intersectData.rayOrig = ray.origin;
        //rtPrintf("made it! %f", hitPt.x);
        // Pass hitPt and normal at hitPt into payload 
        // to calculate payload.radiance in closestHit()

        //payload.hitPoint = hitPt; 
        //    //make_float3(sphere.transform * make_float4(hitPt, 1)); // applying Mp 
        //// transform normal at hitPoint back (M-1)^T
        //payload.hitPointNormal = hitPtNormal;
        //payload.dir = reflectionDir;

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

    //// as explained in: 
    //// https://tavianator.com/2014/ellipsoid_bounding_boxes.html

    //Matrix4x4 S = Matrix4x4::identity(); 
    //float rad_2 = (sphere.radius);
    //S.setRow(0, make_float4(rad_2, .0f, .0f, .0f));
    //S.setRow(1, make_float4(.0f, rad_2, .0f, .0f));
    //S.setRow(2, make_float4(.0f, .0f, rad_2, .0f));
    //S.setRow(3, make_float4(.0f, .0f, .0f, -1.0f));
    ////S[0] = rad_2;
    ////S[5] = rad_2;
    ////S[10] = rad_2;
    ////S[15] = -1.0f;
    //Matrix4x4 M = sphere.transform;
    //////Matrix4x4 Q = M.inverse().transpose() * S * M.inverse();
    ////Matrix4x4 M_T = M.transpose();
    ////Matrix4x4 R = M * S.inverse() * M_T;
    //////Matrix4x4 R = Q.inverse();

    //////if (R == R.transpose())
    //////    rtPrintf("same!");
    //float xMax, xMin, yMax, yMin, zMax, zMin;

    //////xMax = (R[3] + sqrtf(powf(R[3], 2.0f) - (R[15] * R[0]))) / (float)R[15];
    //////xMin = (R[3] - sqrtf(powf(R[3], 2.0f) - (R[15] * R[0]))) / (float)R[15];
    //////yMax = (R[7] + sqrtf(powf(R[7], 2.0f) - (R[15] * R[5]))) / (float)R[15];
    //////yMin = (R[7] - sqrtf(powf(R[7], 2.0f) - (R[15] * R[5]))) / (float)R[15];
    //////zMax = (R[11] + sqrtf(powf(R[11], 2.0f) - (R[10] * R[15]))) / (float)R[15];
    //////zMin = (R[11] - sqrtf(powf(R[11], 2.0f) - (R[10] * R[15]))) / (float)R[15];

    //// we end up only requiring the matrix M in the calculation
    //// after plugging R[i][j] into the plane equations, where 
    //// S is 4x4 mat representing sphere (rad and position), and
    //// R == Q.inverse == M_T.inverse * S * M.inverse. Hence:

    //// eg. xmin or xmax = M[1][4] +/- sqrt( M[1][1]^2 + M[1][2]^2 + M[1][3]^2 ) 
    //xMax = M[3] + sqrtf( powf(M[0], 2.0f) + powf(M[1], 2.0f) + powf(M[2], 2.0f));
    //xMin = M[3] - sqrtf( powf(M[0], 2.0f) + powf(M[1], 2.0f) + powf(M[2], 2.0f));
    //yMax = M[7] + sqrtf( powf(M[4], 2.0f) + powf(M[5], 2.0f) + powf(M[6], 2.0f));
    //yMin = M[7] - sqrtf( powf(M[4], 2.0f) + powf(M[5], 2.0f) + powf(M[6], 2.0f));
    //zMax = M[11] + sqrtf( powf(M[8], 2.0f) + powf(M[9], 2.0f) + powf(M[10], 2.0f));
    //zMin = M[11] - sqrtf( powf(M[8], 2.0f) + powf(M[9], 2.0f) + powf(M[10], 2.0f));

    //result[0] = xMin;
    //result[1] = yMin;
    //result[2] = zMin;
    //result[3] = xMax;
    //result[4] = yMax;
    //result[5] = zMax;
    float x, y, z;    
    x = length(make_float3(sphere.transform.getRow(0)));
    y = length(make_float3(sphere.transform.getRow(1)));
    z = length(make_float3(sphere.transform.getRow(2)));
    result[0] = sphere.transform[3] - x;
    result[1] = sphere.transform[7] - y;
    result[2] = sphere.transform[11] - z;
    result[3] = sphere.transform[3] + x;
    result[4] = sphere.transform[7] + y;
    result[5] = sphere.transform[11] + z;
}