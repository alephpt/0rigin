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
    src/SMFTEngine.cpp
    src/SMFTPipelineFactory.cpp
    src/SMFTBufferManager.cpp
    src/SMFTCompute.cpp
    src/SMFTDescriptorManager.cpp
    src/DiracEvolution.cpp
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

# Stückelberg EM validation test (preserved from Wave1D era)
add_executable(test_stuckelberg_vortex_bfield
    test/test_stuckelberg_vortex_bfield.cpp
)

target_link_libraries(test_stuckelberg_vortex_bfield PRIVATE
    Physics
)

target_include_directories(test_stuckelberg_vortex_bfield PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

set_target_properties(test_stuckelberg_vortex_bfield PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# SMFTEngine + Stückelberg EM integration test
add_executable(test_smft_em_integration
    test/test_smft_em_integration.cpp
    src/SMFTEngine.cpp
    src/SMFTPipelineFactory.cpp
    src/SMFTBufferManager.cpp
    src/SMFTCompute.cpp
    src/SMFTDescriptorManager.cpp
    src/DiracEvolution.cpp
)

target_link_libraries(test_smft_em_integration PRIVATE
    Nova
    Physics
    imgui
    Vulkan::Vulkan
    fftw3f
)

target_include_directories(test_smft_em_integration PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

set_target_properties(test_smft_em_integration PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# DiracEvolution EM coupling test
add_executable(test_dirac_em_coupling
    test/test_dirac_em_coupling.cpp
    src/DiracEvolution.cpp
)

target_link_libraries(test_dirac_em_coupling PRIVATE
    ${FFTW3_LIBRARIES}
)

target_include_directories(test_dirac_em_coupling PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${FFTW3_INCLUDE_DIRS}
)

set_target_properties(test_dirac_em_coupling PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

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
