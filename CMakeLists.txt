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
# MSFT APPLICATION
# ============================================================================

# MSFT sources
set(MSFT_SOURCES
    src/MSFT.cpp
    src/MSFTEngine.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/DiracEvolution.cpp
    main.cpp
)

add_executable(MSFT ${MSFT_SOURCES})

target_include_directories(MSFT PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(MSFT PRIVATE
    Nova
    imgui
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    fftw3f
    pthread
    dl
    m
)

# Set output directory
set_target_properties(MSFT PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# ============================================================================
# TEST EXECUTABLE FOR GPU PIPELINE
# ============================================================================

add_executable(test_msft_gpu
    test/test_msft_gpu.cpp
    src/MSFTEngine.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_msft_gpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(test_msft_gpu PRIVATE
    Nova
    imgui
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    fftw3f
    pthread
    dl
    m
)

set_target_properties(test_msft_gpu PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: Descriptor binding validation
add_executable(test_descriptor_bindings
    test/test_descriptor_bindings.cpp
    src/MSFTEngine.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_descriptor_bindings PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(test_descriptor_bindings PRIVATE
    Nova
    imgui
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    fftw3f
    pthread
    dl
    m
)

set_target_properties(test_descriptor_bindings PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: Stochastic particle evolution
add_executable(test_stochastic_particle
    test/test_stochastic_particle.cpp
    src/MSFTEngine.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_stochastic_particle PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(test_stochastic_particle PRIVATE
    Nova
    imgui
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    fftw3f
    pthread
    dl
    m
)

set_target_properties(test_stochastic_particle PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: Operator splitting validation (infrastructure only)
add_executable(test_operator_splitting
    test/test_operator_splitting.cpp
)

target_include_directories(test_operator_splitting PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
)

set_target_properties(test_operator_splitting PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: Hybrid operator splitting integration test
add_executable(test_hybrid_operator_splitting
    test/test_hybrid_operator_splitting.cpp
    src/MSFTEngine.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_hybrid_operator_splitting PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
)

target_link_libraries(test_hybrid_operator_splitting PRIVATE
    Nova
    imgui
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    fftw3f
    pthread
    dl
    m
)

set_target_properties(test_hybrid_operator_splitting PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: CPU-based operator splitting logic test
add_executable(test_operator_splitting_cpu
    test/test_operator_splitting_cpu.cpp
)

set_target_properties(test_operator_splitting_cpu PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: Full operator splitting simulation (10k steps)
add_executable(test_operator_splitting_full
    test/test_operator_splitting_full.cpp
)

set_target_properties(test_operator_splitting_full PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: CPU-based stochastic validation
add_executable(test_stochastic_cpu
    test/test_stochastic_cpu.cpp
    src/MSFTCommon.cpp
)

target_include_directories(test_stochastic_cpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_stochastic_cpu PRIVATE
    pthread
    m
)

set_target_properties(test_stochastic_cpu PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Comprehensive Dirac-Kuramoto test with full data output
add_executable(test_dirac_stochastic_full
    test/test_dirac_stochastic_full.cpp
    src/MSFTCommon.cpp
)

# Split-operator validation test
add_executable(test_split_operator_validation
    test/test_split_operator_validation.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_dirac_stochastic_full PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_dirac_stochastic_full PRIVATE
    pthread
    m
    stdc++fs
)

set_target_properties(test_dirac_stochastic_full PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

target_include_directories(test_split_operator_validation PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_split_operator_validation PRIVATE
    fftw3f
    pthread
    m
)

set_target_properties(test_split_operator_validation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Phase 2 - Scenario 1: Static Defect Localization (CPU-only)
add_executable(test_scenario1_defect
    test/test_scenario1_defect.cpp
    src/DiracEvolution.cpp
    src/MSFTCommon.cpp
)

target_include_directories(test_scenario1_defect PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_scenario1_defect PRIVATE
    fftw3f
    pthread
    m
)

set_target_properties(test_scenario1_defect PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Ehrenfest theorem validation (critical diagnostic)
add_executable(test_ehrenfest_validation
    test/test_ehrenfest_validation.cpp
    src/DiracEvolution.cpp
    src/MSFTCommon.cpp
)

target_include_directories(test_ehrenfest_validation PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_ehrenfest_validation PRIVATE
    fftw3f
    pthread
    m
)

set_target_properties(test_ehrenfest_validation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Simple warmup test
add_executable(test_simple_warmup
    test/test_simple_warmup.cpp
    src/MSFTEngine.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTCommon.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_simple_warmup PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
)

target_link_libraries(test_simple_warmup PRIVATE
    Nova
    imgui
    Vulkan::Vulkan
    glfw
    pthread
    fftw3f
    m
)

set_target_properties(test_simple_warmup PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Noise sweep test
add_executable(test_noise_sweep
    test/test_noise_sweep.cpp
    src/MSFTEngine.cpp
    src/MSFTBufferManager.cpp
    src/MSFTCompute.cpp
    src/MSFTDescriptorManager.cpp
    src/MSFTPipelineFactory.cpp
    src/MSFTCommon.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_noise_sweep PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
)

target_link_libraries(test_noise_sweep PRIVATE
    Nova
    imgui
    Vulkan::Vulkan
    glfw
    pthread
    fftw3f
    m
)

set_target_properties(test_noise_sweep PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Energy fix test
add_executable(test_energy_fix
    test/test_energy_fix.cpp
    src/MSFTCommon.cpp
    src/DiracEvolution.cpp
)

target_include_directories(test_energy_fix PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_energy_fix PRIVATE
    pthread
    m
    fftw3f
)

set_target_properties(test_energy_fix PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test: Find Critical Sigma
add_executable(test_find_critical_sigma
    test/test_find_critical_sigma.cpp
    src/MSFTCommon.cpp
)

target_include_directories(test_find_critical_sigma PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_find_critical_sigma PRIVATE
    pthread
    m
)

set_target_properties(test_find_critical_sigma PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test: Noise Sweep Corrected
add_executable(test_noise_sweep_corrected
    test/test_noise_sweep_corrected.cpp
    src/MSFTCommon.cpp
)

target_include_directories(test_noise_sweep_corrected PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(test_noise_sweep_corrected PRIVATE
    pthread
    m
)

set_target_properties(test_noise_sweep_corrected PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test: Noise Sweep CPU
add_executable(test_noise_sweep_cpu
    test/test_noise_sweep_cpu.cpp
    src/MSFTCommon.cpp
)

target_include_directories(test_noise_sweep_cpu PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(test_noise_sweep_cpu PRIVATE
    ${Vulkan_LIBRARIES}
    pthread
    m
)

set_target_properties(test_noise_sweep_cpu PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# ============================================================================
# SHADERS
# ============================================================================

# Copy shaders to build directory
add_custom_command(TARGET MSFT POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/shaders
    ${CMAKE_BINARY_DIR}/shaders
)

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS "===============================================")
message(STATUS " 0rigin: MSFT+Dirac GPU Simulation")
message(STATUS "===============================================")
message(STATUS "")
message(STATUS "  Vulkan SDK: ${Vulkan_INCLUDE_DIRS}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "  Targets:")
message(STATUS "    • Nova (static library)")
message(STATUS "    • imgui (static library)")
message(STATUS "    • MSFT (executable)")
message(STATUS "")
message(STATUS "===============================================")