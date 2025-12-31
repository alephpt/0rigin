cmake_minimum_required(VERSION 3.20)
project(0rigin-wave1d LANGUAGES CXX)

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
# PHYSICS LIBRARY
# ============================================================================

add_library(Physics STATIC
    src/physics/ProcaEM.cpp
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
# SUMMARY
# ============================================================================

message(STATUS "")
message(STATUS "===============================================")
message(STATUS " 0rigin: Wave1D Gauge-Covariant EM")
message(STATUS "===============================================")
message(STATUS "")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "")
message(STATUS "  Targets:")
message(STATUS "    • Physics (static library)")
message(STATUS "")
message(STATUS "===============================================")
