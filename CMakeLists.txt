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
    src/physics/ProcaVulkanPipeline.cpp
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
# SMFTCORE LIBRARY
# ============================================================================

add_library(SMFTCore STATIC
    src/SMFTCore.cpp
)

target_include_directories(SMFTCore PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(SMFTCore PUBLIC
    Physics
)

target_compile_options(SMFTCore PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

# ============================================================================
# SMFT EXECUTABLE (Minimal Demo)
# ============================================================================

add_executable(smft
    src/main_smft.cpp
)

target_include_directories(smft PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(smft PRIVATE
    SMFTCore
    Physics
    pthread
    m
)

target_compile_options(smft PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)

set_target_properties(smft PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
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
message(STATUS "    • SMFTCore (static library)")
message(STATUS "    • smft (executable)")
message(STATUS "")
message(STATUS "===============================================")
