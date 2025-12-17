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
    pthread
    dl
    m
)

set_target_properties(test_stochastic_particle PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Test executable: CPU-based stochastic validation
add_executable(test_stochastic_cpu
    test/test_stochastic_cpu.cpp
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