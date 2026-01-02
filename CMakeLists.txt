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
# SMFT APPLICATION (GPU-accelerated + Stückelberg EM)
# ============================================================================

# SMFT sources - unified executable with interactive and test modes
set(SMFT_SOURCES
    src/SMFTCore.cpp
    src/SMFTCore3D.cpp
    src/SMFTEngine.cpp
    src/SMFTEngine3D.cpp
    src/SMFTPipelineFactory.cpp
    src/SMFTBufferManager.cpp
    src/SMFTCompute.cpp
    src/SMFTDescriptorManager.cpp
    src/DiracEvolution.cpp
    src/GeodesicIntegrator.cpp
    src/output/OutputManager.cpp
    src/simulations/TestConfig.cpp
    src/simulations/SMFTTestRunner.cpp
    src/simulations/ObservableComputer.cpp
    main.cpp
)

add_executable(SMFT ${SMFT_SOURCES})

target_include_directories(SMFT PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(SMFT PRIVATE
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

# Set output directory and binary name
set_target_properties(SMFT PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    OUTPUT_NAME "smft"
)

# ============================================================================
# TEST EXECUTABLES
# ============================================================================

# Geodesic Equation Verification Test
add_executable(test_geodesic_verification
    test/test_geodesic_verification.cpp
    src/GeodesicIntegrator.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_geodesic_verification PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
)

target_link_libraries(test_geodesic_verification PRIVATE
    fftw3f
    pthread
    m
)

set_target_properties(test_geodesic_verification PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# A2: Weak Field Limit Validation Test
add_executable(test_weak_field_limit
    test/test_weak_field_limit.cpp
    src/simulations/SMFTTestRunner.cpp
    src/simulations/TestConfig.cpp
    src/simulations/ObservableComputer.cpp
    src/SMFTCore.cpp
    src/SMFTEngine.cpp
    src/SMFTPipelineFactory.cpp
    src/SMFTBufferManager.cpp
    src/SMFTCompute.cpp
    src/SMFTDescriptorManager.cpp
    src/DiracEvolution.cpp
    src/GeodesicIntegrator.cpp
    src/output/OutputManager.cpp
)

target_include_directories(test_weak_field_limit PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(test_weak_field_limit PRIVATE
    Nova
    Physics
    fftw3f
    yaml-cpp
    pthread
    m
    ${Vulkan_LIBRARIES}
)

set_target_properties(test_weak_field_limit PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# SMFTCore3D Basic Test
add_executable(test_smftcore3d_basic
    test/test_smftcore3d_basic.cpp
    src/SMFTCore3D.cpp
)

target_include_directories(test_smftcore3d_basic PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(test_smftcore3d_basic PRIVATE
    ${Vulkan_LIBRARIES}
    pthread
    m
)

set_target_properties(test_smftcore3d_basic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# SMFTCore3D GPU Test (Simplified)
add_executable(test_smftcore3d_gpu
    test/test_smftcore3d_gpu_simple.cpp
    src/SMFTCore3D.cpp
)

target_include_directories(test_smftcore3d_gpu PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(test_smftcore3d_gpu PRIVATE
    ${Vulkan_LIBRARIES}
    pthread
    m
)

set_target_properties(test_smftcore3d_gpu PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# SMFT 3D CPU-Only Test (no GPU dependency)
add_executable(test_smft3d_cpu_only
    test/test_smft3d_cpu_only.cpp
    src/SMFTCore3D.cpp
)

target_include_directories(test_smft3d_cpu_only PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(test_smft3d_cpu_only PRIVATE
    ${Vulkan_LIBRARIES}
    pthread
    m
)

set_target_properties(test_smft3d_cpu_only PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# ============================================================================
# LEGACY TEST BINARIES REMOVED
# ============================================================================
# All tests now run through unified smft binary via YAML configs:
#   ./smft --test config/stuckelberg_integration_test.yaml
#
# Removed binaries:
# - test_stuckelberg_vortex_bfield → config/stuckelberg_integration_test.yaml
# - test_smft_em_integration → integrated into smft
# - test_dirac_em_coupling → integrated into smft
#
# Single source of truth: ./smft

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS "===============================================")
message(STATUS " 0rigin: GPU-Accelerated SMFT + Stückelberg EM")
message(STATUS "===============================================")
message(STATUS "")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "  Libraries:")
message(STATUS "    • Physics (Stückelberg EM)")
message(STATUS "    • Nova (Vulkan compute engine)")
message(STATUS "    • ImGui (UI framework)")
message(STATUS "")
message(STATUS "  Targets:")
message(STATUS "    • smft (main executable)")
message(STATUS "       - GPU-accelerated Kuramoto dynamics")
message(STATUS "       - CPU Dirac evolution (split-operator)")
message(STATUS "       - Stückelberg EM coupling (ready to integrate)")
message(STATUS "    • test_stuckelberg_vortex_bfield (EM validation)")
message(STATUS "    • test_dirac_em_coupling (Dirac EM field coupling)")
message(STATUS "")
message(STATUS "===============================================")
