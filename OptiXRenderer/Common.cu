#include <optix.h>
#include <optix_device.h>

#include "Payloads.h"

using namespace optix;

rtDeclareVariable(Payload, payload, rtPayload, );
rtDeclareVariable(float3, backgroundColor, , );

RT_PROGRAM void miss()
{
    // Set the result to be the background color if miss
    // TODO: change the color to backgroundColor
    //backgroundColor = make_float3(0, 0, 1);
    //payload.radiance = make_float3(1, 0, 0);
    payload.radiance = backgroundColor; 
    payload.done = true;
    //rtPrintf("miss!");
}

RT_PROGRAM void exception()
{
    // Print any exception for debugging
    const unsigned int code = rtGetExceptionCode();
    rtPrintExceptionDetails();
}

rtDeclareVariable(ShadowPayload, shadowPayload, rtPayload, );
rtDeclareVariable(float1, t, rtIntersectionDistance, );

RT_PROGRAM void anyHit()
{
    //rtPrintf("anyHit!");
    shadowPayload.isVisible = false;
    rtTerminateRay();
}