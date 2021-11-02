set(CPP_FILES 
    main.cpp
    Renderer.cpp
    SceneLoader.cpp
    utils.cpp
    lodepng.cpp
)

set(H_FILES
    main.h
    Renderer.h 
    SceneLoader.h 
    utils.h 
    lodepng.h
    Scene.h
    Light.h
    Geometries.h
    Payloads.h
<<<<<<< HEAD
    Config.h 
=======
>>>>>>> f6caa29044f36949ca28b0fb6334733497a0497d
    random.h
)

set(CU_FILES
    PinholeCamera.cu 
    RayTracer.cu 
    Triangle.cu
    Sphere.cu 
    Common.cu
) 

include_directories(${GLUT_INCLUDE_DIR})

# See top level CMakeLists.txt file for documentation of OPTIX_add_sample_executable.
OPTIX_add_sample_executable(OptiXRenderer
  ${CPP_FILES}
  ${H_FILES}
  ${CU_FILES} 
)