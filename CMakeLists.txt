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

# Python3 for matplotlib-cpp
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
if(NOT Python3_FOUND)
    message(FATAL_ERROR "Python3 with NumPy not found")
endif()
message(STATUS "Python3: ${Python3_VERSION}")
message(STATUS "Python3 Include: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 Libraries: ${Python3_LIBRARIES}")

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
# SMFT APPLICATION
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
    src/KleinGordonEvolution.cpp
    src/SMFTCommon.cpp
    src/output/OutputManager.cpp
    src/output/PlottingManager.cpp
    src/simulations/TestConfig.cpp
    src/simulations/SMFTTestRunner.cpp
    src/simulations/ObservableComputer.cpp
    src/simulations/EMObservables.cpp
    src/simulations/InitialConditions.cpp
    src/validation/ValidationCriteria.cpp
    src/validation/GlobalValidator.cpp
    src/validation/ScenarioValidator.cpp
    src/validation/PhaseTransitionAnalyzer.cpp
    src/validation/EnergyBudget.cpp
    src/physics/LorentzTransform.cpp
    src/physics/EMFieldComputer.cpp
    src/analysis/PowerLawFitter.cpp
    src/analysis/GeometryAnalyzer.cpp
    src/analysis/TrajectoryComparator.cpp
    src/analysis/DualSimulationComparator.cpp
    src/analysis/ForceDecomposer.cpp
    main.cpp
)

add_executable(SMFT ${SMFT_SOURCES})

target_include_directories(SMFT PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/lib
    ${CMAKE_SOURCE_DIR}/lib/Nova
    ${CMAKE_SOURCE_DIR}/lib/matplotlib-cpp
    ${Vulkan_INCLUDE_DIRS}
    ${SDL2_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${Python3_NumPy_INCLUDE_DIRS}
)

target_link_libraries(SMFT PRIVATE
    Nova
    imgui
    ${Vulkan_LIBRARIES}
    ${SDL2_LIBRARIES}
    ${Python3_LIBRARIES}
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
# TEST FRAMEWORK - YAML-DRIVEN TESTS ONLY
# ============================================================================
# All tests MUST use the YAML-driven test framework via SMFTTestRunner.
# Run tests using: ./bin/smft --test config/<test_name>.yaml
#
# Single Source of Truth: SMFTEngine (GPU shaders + Dirac CPU evolution)
# - NO standalone test binaries (violates consistency)
# - NO MSFTCommon CPU reimplementation (wrong physics)
# - ALL tests validate the same SMFTEngine implementation
#
# Example configs:
#   ./bin/smft --test config/timesync_validation.yaml
#   ./bin/smft --test config/traveling_wave_validation.yaml
#   ./bin/smft --test config/quick_validation.yaml
# ============================================================================

# ============================================================================
# SHADERS
# ============================================================================

# Copy shaders to build directory
add_custom_command(TARGET SMFT POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/shaders
    ${CMAKE_BINARY_DIR}/shaders
)

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS " 0rigin: SMFT+Dirac GPU Simulation")
message(STATUS "===============================================")
message(STATUS "")
message(STATUS "  Vulkan SDK: ${Vulkan_INCLUDE_DIRS}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "  Targets:")
message(STATUS "    • Nova (static library)")
message(STATUS "    • imgui (static library)")
message(STATUS "    • SMFT (executable)")
message(STATUS "")
message(STATUS "===============================================")