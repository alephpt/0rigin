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
    src/Maxwell3D.cpp
    src/Dirac3D.cpp
    src/TRDEngine.cpp
    src/TRDEngine3D.cpp
    src/TRDPipelineFactory.cpp
    src/TRDBufferManager.cpp
    src/TRDCompute.cpp
    src/TRDDescriptorManager.cpp
    src/DiracEvolution.cpp
    # src/FieldEvolution.cpp  # Commented out - file doesn't exist yet
    src/GeodesicIntegrator.cpp
    src/output/OutputManager.cpp
    src/simulations/TestConfig.cpp
    src/simulations/TRDTestRunner.cpp
    src/simulations/ObservableComputer.cpp
    # 3D particle physics tests (integrated into single executable)
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
    # B3-B6 Standard Model tests
    test/test_three_generations.cpp
    test/test_electroweak.cpp
    test/test_strong_force.cpp
    test/test_higgs_connection.cpp
    # D2 Hardware experimental tests
    test/test_josephson_junction.cpp
    # D3 Gravitational wave tests
    test/test_binary_merger.cpp
    # D4 Particle scattering test
    test/test_particle_scattering.cpp
    # H3 Spin-Magnetism test
    test/test_spin_magnetism.cpp
    # H3 Magnetic Dynamo test
    test/test_magnetic_dynamo.cpp
    # H1 Knot Topology test
    test/test_knot_topology.cpp
    test/test_knot_stability.cpp
    # H2 Solar System Orbital Dynamics test
    test/test_solar_system.cpp
    # F2 Multi-Scale Validation test
    test/test_multiscale.cpp
    # F3 Finite Temperature Effects test
    test/test_finite_temperature.cpp
    # F4 Quantum Fluctuation Incorporation test
    test/test_quantum_fluctuations.cpp
    # F5 High-Performance Computing Scaling test
    test/test_hpc_scaling.cpp
    # B2 Fine Structure Constant test
    test/test_fine_structure_constant.cpp
    # D5 Atomic Physics and Precision Spectroscopy test
    test/test_atomic_physics.cpp
    # D3 Astrophysical Observations test
    test/test_astrophysical_observations.cpp
    # D4 LHC Predictions test
    test/test_lhc_predictions.cpp
    # G4 Navier-Stokes Fluid Dynamics test
    test/test_navier_stokes.cpp
    # E1/E3/E5 Wave 1 Mathematical Rigor tests
    test/test_renormalizability.cpp
    test/test_causality.cpp
    # D1 Solar System test
    test/test_solar_system.cpp
    # H3 Magnetic Dynamo test
    test/test_magnetic_dynamo.cpp
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
# TEST EXECUTABLES
# ============================================================================
# IMPORTANT: CLAUDE.md ABSOLUTE RULE - Only ./trd --test <config.yaml> allowed
# Most tests have been migrated to unified TRD executable
# Only unit tests for core infrastructure remain as standalone binaries
# ============================================================================

# D1: Experimental Predictions Test (Falsifiability)
# MIGRATED to ./trd --test config/experimental_predictions.yaml
# add_executable(test_experimental_predictions
#     test/test_experimental_predictions.cpp
# )
# target_include_directories(test_experimental_predictions PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
# )
# target_link_libraries(test_experimental_predictions PRIVATE
#     m
# )
# set_target_properties(test_experimental_predictions PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# Geodesic Equation Verification Test
# MIGRATED to ./trd --test config/geodesic_3d.yaml
# add_executable(test_geodesic_verification
#     test/test_geodesic_verification.cpp
#     src/GeodesicIntegrator.cpp
#     src/DiracEvolution.cpp
# )
# target_include_directories(test_geodesic_verification PRIVATE
#     ${CMAKE_SOURCE_DIR}/src
#     ${CMAKE_SOURCE_DIR}/include
#     ${CMAKE_SOURCE_DIR}/lib
# )
# target_link_libraries(test_geodesic_verification PRIVATE
#     fftw3f
#     pthread
#     m
# )
# set_target_properties(test_geodesic_verification PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# E1: Renormalizability Analysis Test (Mathematical Rigor)
# REMOVED - Integrated into unified trd executable
# Use: ./trd --test config/renormalizability.yaml
# add_executable(test_renormalizability
#     test/test_renormalizability.cpp
# )
# target_include_directories(test_renormalizability PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
# )
# target_link_libraries(test_renormalizability PRIVATE
#     yaml-cpp
#     m
# )
# set_target_properties(test_renormalizability PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# E3: Causality Validation Test
# REMOVED - Integrated into unified trd executable
# Use: ./trd --test config/causality.yaml
# add_executable(test_causality
#     test/test_causality.cpp
#     src/TRDCore3D.cpp
# )
# target_include_directories(test_causality PRIVATE
#     ${CMAKE_SOURCE_DIR}/src
#     ${CMAKE_SOURCE_DIR}/include
#     ${Vulkan_INCLUDE_DIRS}
# )
# target_link_libraries(test_causality PRIVATE
#     ${Vulkan_LIBRARIES}
#     pthread
#     m
# )
# set_target_properties(test_causality PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# A1: Einstein Field Equations Test
# MIGRATED to ./trd --test (integrated into main.cpp)
# add_executable(test_einstein_field_equations
#     test/test_einstein_field_equations.cpp
#     src/TRDCore3D.cpp
#     src/Maxwell3D.cpp
# )
# target_include_directories(test_einstein_field_equations PRIVATE
#     ${CMAKE_SOURCE_DIR}/src
#     ${CMAKE_SOURCE_DIR}/include
# )
# target_link_libraries(test_einstein_field_equations PRIVATE
#     ${CMAKE_THREAD_LIBS_INIT}
#     m
# )
# set_target_properties(test_einstein_field_equations PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# A2: Weak Field Limit Validation Test
# MIGRATED to ./trd --test config/weak_field_3d.yaml
# add_executable(test_weak_field_limit
#     test/test_weak_field_limit.cpp
#     src/simulations/TRDTestRunner.cpp
#     src/simulations/TestConfig.cpp
#     src/simulations/ObservableComputer.cpp
#     src/TRDCore.cpp
#     src/TRDEngine.cpp
#     src/TRDPipelineFactory.cpp
#     src/TRDBufferManager.cpp
#     src/TRDCompute.cpp
#     src/TRDDescriptorManager.cpp
#     src/DiracEvolution.cpp
#     src/GeodesicIntegrator.cpp
#     src/output/OutputManager.cpp
# )
# target_include_directories(test_weak_field_limit PRIVATE
#     ${CMAKE_SOURCE_DIR}/src
#     ${CMAKE_SOURCE_DIR}/include
#     ${CMAKE_SOURCE_DIR}/lib
#     ${CMAKE_SOURCE_DIR}/lib/Nova
#     ${Vulkan_INCLUDE_DIRS}
#     ${SDL2_INCLUDE_DIRS}
# )
# target_link_libraries(test_weak_field_limit PRIVATE
#     Nova
#     Physics
#     fftw3f
#     yaml-cpp
#     pthread
#     m
#     ${Vulkan_LIBRARIES}
# )
# set_target_properties(test_weak_field_limit PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# TRDCore3D Basic Test
add_executable(test_trdcore3d_basic
    test/test_trdcore3d_basic.cpp
    src/TRDCore3D.cpp
)

target_include_directories(test_trdcore3d_basic PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(test_trdcore3d_basic PRIVATE
    ${Vulkan_LIBRARIES}
    pthread
    m
)

set_target_properties(test_trdcore3d_basic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# TRDCore3D GPU Test (Simplified)
# MIGRATED to ./trd --test (GPU functionality tested in main executable)
# add_executable(test_trdcore3d_gpu
#     test/test_trdcore3d_gpu_simple.cpp
#     src/TRDCore3D.cpp
# )
# target_include_directories(test_trdcore3d_gpu PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
#     ${CMAKE_SOURCE_DIR}/src
#     ${CMAKE_SOURCE_DIR}/lib
#     ${CMAKE_SOURCE_DIR}/lib/Nova
#     ${Vulkan_INCLUDE_DIRS}
# )
# target_link_libraries(test_trdcore3d_gpu PRIVATE
#     ${Vulkan_LIBRARIES}
#     pthread
#     m
# )
# set_target_properties(test_trdcore3d_gpu PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# TRDCore3D Symplectic Integration Test
add_executable(test_trdcore3d_symplectic
    test/test_trdcore3d_symplectic.cpp
    src/TRDCore3D.cpp
)

target_include_directories(test_trdcore3d_symplectic PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/lib
    ${Vulkan_INCLUDE_DIRS}
)

target_link_libraries(test_trdcore3d_symplectic PRIVATE
    ${Vulkan_LIBRARIES}
    pthread
    m
)

set_target_properties(test_trdcore3d_symplectic PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# TRD 3D CPU-Only Test (no GPU dependency)
# MIGRATED to ./trd --test (CPU mode available in main executable)
# add_executable(test_trd3d_cpu_only
#     test/test_trd3d_cpu_only.cpp
#     src/TRDCore3D.cpp
# )
# target_include_directories(test_trd3d_cpu_only PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
#     ${Vulkan_INCLUDE_DIRS}
# )
# target_link_libraries(test_trd3d_cpu_only PRIVATE
#     ${Vulkan_LIBRARIES}
#     pthread
#     m
# )
# set_target_properties(test_trd3d_cpu_only PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# Maxwell3D Wave Propagation Test (Week 5-6)
# MIGRATED to ./trd --test (Maxwell3D integrated into main)
# add_executable(test_maxwell3d_wave
#     test/test_maxwell3d_wave.cpp
#     src/Maxwell3D.cpp
# )
# target_include_directories(test_maxwell3d_wave PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
# )
# target_link_libraries(test_maxwell3d_wave PRIVATE
#     pthread
#     m
# )
# set_target_properties(test_maxwell3d_wave PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# G1: EM Wave Propagation 3D Test (comprehensive wave velocity validation)
# MIGRATED to ./trd --test (EM wave tests integrated)
# add_executable(test_em_wave_propagation_3d
#     test/test_em_wave_propagation_3d.cpp
#     src/Maxwell3D.cpp
# )
# target_include_directories(test_em_wave_propagation_3d PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
# )
# target_link_libraries(test_em_wave_propagation_3d PRIVATE
#     pthread
#     m
# )
# set_target_properties(test_em_wave_propagation_3d PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# ============================================================================
# 2D→3D VALIDATION MIGRATION TEST SUITE
# ============================================================================
#
# CONSOLIDATED: All 3D particle physics tests now run through single ./trd executable
#   Usage: ./trd --test config/<test_name>_3d.yaml
#
# Integrated tests:
#   - Lorentz Force 3D        → ./trd --test config/lorentz_force_3d.yaml
#   - Stückelberg Vortex 3D   → ./trd --test config/stuckelberg_vortex_3d.yaml
#   - Geodesic 3D             → ./trd --test config/geodesic_3d.yaml
#   - Weak Field 3D           → ./trd --test config/weak_field_3d.yaml
#   - Three-Body EM 3D        → ./trd --test config/three_body_em_3d.yaml
#   - EM-Gravity Coupling 3D  → ./trd --test config/em_gravity_coupling_3d.yaml
#
# REMOVED standalone test executables (architecture violation - single executable only)

# Dirac3D Free Particle Test (Week 7-8)
# MIGRATED to ./trd --test (Dirac3D integrated into main)
# add_executable(test_dirac3d_free
#     test/test_dirac3d_free.cpp
#     src/Dirac3D.cpp
# )
# target_include_directories(test_dirac3d_free PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
# )
# target_link_libraries(test_dirac3d_free PRIVATE
#     fftw3f
#     pthread
#     m
# )
# set_target_properties(test_dirac3d_free PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# Full 3D TRD Integration Test (Week 9-10)
# MIGRATED to ./trd --test (full 3D integration in main executable)
# add_executable(test_3d_full_trd
#     test/test_3d_full_trd.cpp
#     src/TRDCore3D.cpp
#     src/Maxwell3D.cpp
#     src/Dirac3D.cpp
# )
# target_include_directories(test_3d_full_trd PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
#     ${Vulkan_INCLUDE_DIRS}
# )
# target_link_libraries(test_3d_full_trd PRIVATE
#     ${Vulkan_LIBRARIES}
#     fftw3f
#     pthread
#     m
# )
# set_target_properties(test_3d_full_trd PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# FieldEvolution - Symplectic Integrators Test
# Commented out - files don't exist yet
# add_executable(test_field_evolution
#     test/test_field_evolution.cpp
#     src/FieldEvolution.cpp
# )
#
# target_include_directories(test_field_evolution PRIVATE
#     ${CMAKE_SOURCE_DIR}/include
# )
#
# target_link_libraries(test_field_evolution PRIVATE
#     pthread
#     m
# )
#
# set_target_properties(test_field_evolution PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
# )

# Field Initializers Validation - Verify TRDFieldInitializers.h
add_executable(test_field_initializers_validation
    test/test_field_initializers_validation.cpp
)

target_include_directories(test_field_initializers_validation PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_field_initializers_validation PRIVATE
    m
)

set_target_properties(test_field_initializers_validation PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# ============================================================================
# LEGACY TEST BINARIES REMOVED
# ============================================================================
# All tests now run through unified trd binary via YAML configs:
#   ./trd --test config/stuckelberg_integration_test.yaml
#
# Removed binaries:
# - test_stuckelberg_vortex_bfield → config/stuckelberg_integration_test.yaml
# - test_trd_em_integration → integrated into trd
# - test_dirac_em_coupling → integrated into trd
#
# Single source of truth: ./trd

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS "===============================================")
message(STATUS " 0rigin: GPU-Accelerated TRD + Stückelberg EM")
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
message(STATUS "    • trd (main executable)")
message(STATUS "       - GPU-accelerated Kuramoto dynamics")
message(STATUS "       - CPU Dirac evolution (split-operator)")
message(STATUS "       - Stückelberg EM coupling (ready to integrate)")
message(STATUS "    • test_stuckelberg_vortex_bfield (EM validation)")
message(STATUS "    • test_dirac_em_coupling (Dirac EM field coupling)")
message(STATUS "")
message(STATUS "===============================================")
