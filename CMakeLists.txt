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

# ============================================================================
# SUBDIRECTORIES
# ============================================================================

# Nova library
add_subdirectory(lib/Nova)

# Nova SMFT library
add_subdirectory(lib/Nova/Core/smft)

# Examples
add_subdirectory(examples/vulkan)

# ============================================================================
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS "===============================================")
message(STATUS " 0rigin: SMFT+Dirac GPU Simulation")
message(STATUS "===============================================")
message(STATUS "")
message(STATUS "  Vulkan SDK: ${Vulkan_INCLUDE_DIRS}")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "  Targets:")
message(STATUS "    • NovaSMFT (static library)")
message(STATUS "    • smft_test (executable)")
message(STATUS "")
message(STATUS "===============================================")
