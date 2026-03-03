cmake_minimum_required(VERSION 3.20)
project(0rigin LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ============================================================================
# DEPENDENCIES
# ============================================================================

find_package(Vulkan REQUIRED)

if(NOT Vulkan_FOUND)
    message(FATAL_ERROR "Vulkan SDK not found. Install from https://vulkan.lunarg.com/")
endif()

message(STATUS "Vulkan SDK: ${Vulkan_INCLUDE_DIRS}")

# SDL2
find_package(SDL2 REQUIRED)
if(NOT SDL2_FOUND)
    message(FATAL_ERROR "SDL2 not found")
endif()

# FFTW3
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3f)
if(NOT FFTW3_FOUND)
    message(FATAL_ERROR "FFTW3 not found")
endif()

# OpenMP (optional, for HPC scaling tests)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
else()
    message(STATUS "OpenMP not found - HPC scaling tests will run single-threaded")
endif()

# ============================================================================
# PHYSICS LIBRARY (Stückelberg EM)
# ============================================================================

add_library(Physics STATIC
    src/physics/StuckelbergEM.cpp
)

target_include_directories(Physics PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_compile_options(Physics PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

# ============================================================================
# NOVA LIBRARY
# ============================================================================

# Collect Nova source files
file(GLOB_RECURSE NOVA_SOURCES
    lib/Nova/*.cpp
    lib/Nova/Core/*.cpp
    lib/Nova/Core/modules/*.cpp
    lib/Nova/Core/modules/pipeline/*.cpp
    lib/Nova/Core/modules/camera/*.cpp
    lib/Nova/Core/modules/atomic/*.cpp
    lib/Nova/Core/components/*.cpp
)

# Remove example and test files from Nova sources
list(FILTER NOVA_SOURCES EXCLUDE REGEX ".*Example.*")
list(FILTER NOVA_SOURCES EXCLUDE REGEX ".*test.*")
list(FILTER NOVA_SOURCES EXCLUDE REGEX ".*extern/fastgltf.*")

# Nova library
add_library(Nova STATIC ${NOVA_SOURCES})

target_include_directories(Nova PUBLIC
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${CMAKE_SOURCE_DIR}/lib/Nova/Core
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(Nova PUBLIC
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
)

# ============================================================================
# IMGUI LIBRARY
# ============================================================================

file(GLOB IMGUI_SOURCES
    lib/imgui/*.cpp
    lib/imgui/backends/imgui_impl_sdl2.cpp
    lib/imgui/backends/imgui_impl_vulkan.cpp
)

add_library(imgui STATIC ${IMGUI_SOURCES})

target_include_directories(imgui PUBLIC
    ${CMAKE_SOURCE_DIR}/lib/imgui
    ${CMAKE_SOURCE_DIR}/lib/imgui/backends
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

# ============================================================================
# TRD APPLICATION (GPU-accelerated + Stückelberg EM)
# ============================================================================

# TRD sources - unified executable with interactive and test modes
set(TRD_SOURCES
    src/TRDCore.cpp
    src/TRDCore3D.cpp
    src/ConservativeSolver.cpp
    src/Maxwell3D.cpp
    src/Dirac3D.cpp
    src/TRDEngine.cpp
    src/TRDEngine3D.cpp
    src/TRDPipelineFactory.cpp
    src/TRDBufferManager.cpp
    src/TRDCompute.cpp
    src/TRDDescriptorManager.cpp
    src/DiracEvolution.cpp
    src/GeodesicIntegrator.cpp
    src/output/OutputManager.cpp
    src/simulations/TestConfig.cpp
    src/simulations/TRDTestRunner.cpp
    src/simulations/ObservableComputer.cpp
    test/test_lorentz_force_3d.cpp
    test/test_three_body_em_3d.cpp
    test/test_stuckelberg_vortex_3d.cpp
    test/test_geodesic_3d.cpp
    test/test_weak_field_3d.cpp
    test/test_time_dilation_3d.cpp
    test/test_gravitational_waves.cpp
    test/test_em_gravity_coupling_3d.cpp
    test/test_einstein_field_equations.cpp
    test/test_light_deflection_3d.cpp
    test/test_particle_spectrum_unified.cpp
    test/test_cosmological_constant.cpp
    test/test_friedmann_equations.cpp
    test/test_dark_matter.cpp
    test/test_dark_energy.cpp
    test/test_inflation.cpp
    test/test_unitarity.cpp
    test/test_scale_invariance.cpp
    test/test_symmetry_analysis.cpp
    test/test_three_generations.cpp
    test/test_electroweak.cpp
    test/test_strong_force.cpp
    test/test_higgs_connection.cpp
    test/test_laboratory_scale.cpp
    test/test_josephson_junction.cpp
    test/test_binary_merger.cpp
    test/test_particle_scattering.cpp
    test/test_spin_magnetism.cpp
    test/test_magnetic_dynamo.cpp
    test/test_knot_topology.cpp
    test/test_knot_stability.cpp
    test/test_solar_system.cpp
    test/test_multiscale.cpp
    test/test_finite_temperature.cpp
    test/test_quantum_fluctuations.cpp
    test/test_hpc_scaling.cpp
    test/test_fine_structure_constant.cpp
    test/test_atomic_physics.cpp
    test/test_astrophysical_observations.cpp
    test/test_lhc_predictions.cpp
    test/test_navier_stokes.cpp
    test/test_renormalizability.cpp
    test/test_causality.cpp
    main.cpp
)

add_executable(TRD ${TRD_SOURCES})

target_compile_definitions(TRD PRIVATE TRD_MAIN_EXECUTABLE)

target_include_directories(TRD PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(TRD PRIVATE
    Nova
    imgui
    Physics
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    fftw3f
    yaml-cpp
    pthread
    dl
    m
)

# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(TRD PRIVATE OpenMP::OpenMP_CXX)
endif()

# Set output directory and binary name
set_target_properties(TRD PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    OUTPUT_NAME "trd"
)

# ============================================================================
# SHADER COMPILATION (SPIR-V from GLSL)
# ============================================================================

find_program(GLSLC glslc)
if(NOT GLSLC)
    message(WARNING "glslc not found. Shader compilation will be skipped.")
    message(WARNING "Install glslc: sudo apt-get install glslang-tools")
else()
    message(STATUS "Found glslc: ${GLSLC}")

    # Find all GLSL compute shaders
    file(GLOB_RECURSE SHADER_SOURCES
        "${CMAKE_SOURCE_DIR}/shaders/smft/*.comp"
    )

    set(COMPILED_SHADERS "")

    # Compile each shader to SPIR-V
    foreach(SHADER_SOURCE ${SHADER_SOURCES})
        get_filename_component(SHADER_NAME ${SHADER_SOURCE} NAME)
        get_filename_component(SHADER_DIR ${SHADER_SOURCE} DIRECTORY)
        set(SHADER_OUTPUT "${SHADER_DIR}/${SHADER_NAME}.spv")

        add_custom_command(
            OUTPUT ${SHADER_OUTPUT}
            COMMAND ${GLSLC} ${SHADER_SOURCE} -o ${SHADER_OUTPUT}
            DEPENDS ${SHADER_SOURCE}
            COMMENT "Compiling shader: ${SHADER_NAME} -> ${SHADER_NAME}.spv"
            VERBATIM
        )

        list(APPEND COMPILED_SHADERS ${SHADER_OUTPUT})
    endforeach()

    # Create shader compilation target
    if(COMPILED_SHADERS)
        add_custom_target(CompileShaders ALL DEPENDS ${COMPILED_SHADERS})

        # Ensure shaders compile before TRD binary
        add_dependencies(TRD CompileShaders)

        list(LENGTH COMPILED_SHADERS SHADER_COUNT)
        message(STATUS "Shader compilation configured: ${SHADER_COUNT} shaders from ${CMAKE_SOURCE_DIR}/shaders/")
    endif()
endif()

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS "===============================================")
message(STATUS " 0rigin: GPU-Accelerated TRD + Stuckelberg EM")
message(STATUS "===============================================")
message(STATUS "")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "  Libraries:")
message(STATUS "    - Physics (Stuckelberg EM)")
message(STATUS "    - Nova (Vulkan compute engine)")
message(STATUS "    - ImGui (UI framework)")
message(STATUS "")
message(STATUS "  Targets:")
message(STATUS "    - trd (unified executable)")
message(STATUS "      Usage: ./trd --test config/<test>.yaml")
message(STATUS "")
message(STATUS "===============================================")
