#ifndef _CAMERA_H_
#define _CAMERA_H_

#include <optixu/optixu_math_namespace.h>
#include <optixu/optixu_matrix_namespace.h>

struct Camera {

	optix::float3 pos; // lookFrom
	optix::float3 lookAt; 
	optix::float3 up;
	float fov;

	// TODO: find w, v and u using cross products
	optix::float3 w; 
	optix::float3 u; 
	optix::float3 v; 

	Camera(optix::float3 look_from, 
		optix::float3 look_at, optix::float3 up_vec, float field_of_view) {

		pos = look_from;
		lookAt = look_at;
		up = up_vec;
		fov = field_of_view;

		w = optix::normalize(pos - lookAt);
		optix::float3 b = up;
		u = optix::normalize(optix::cross(b, w));
		v = optix::cross(w, u); // since w and u are normalize, v is normalized
	}

	Camera() {
		pos = optix::make_float3(.0f);
		lookAt = optix::make_float3(.0f);
		up = optix::make_float3(.0f);
		fov = .0f;

		w = optix::make_float3(.0f);
		u = optix::make_float3(.0f);
		v = optix::make_float3(.0f);
	}

};

#endif