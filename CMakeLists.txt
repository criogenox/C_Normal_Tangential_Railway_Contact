cmake_minimum_required(VERSION 4.2)
project(ContactSuite LANGUAGES CXX)

# Precompiled OpenMP library / Use FetchContent module instead
set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Release as default
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            Debug Release RelWithDebInfo MinSizeRel)
endif ()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion)
endif ()

option(BUILD_EQ_CONTACT "Build Eq-Elastic-Norm module" ON)
option(BUILD_TABCON "Build tabcon module" ON)
option(BUILD_FASTSIM "Build FastSim module" ON)

if (BUILD_EQ_CONTACT)
    message(STATUS "Building Equivalent Elastic Normal module...")
    add_subdirectory(normal_eq_elastic)
endif ()

if (BUILD_TABCON)
    message(STATUS "Building tabcon module...")
    add_subdirectory(tang_spin_tabcon)
endif ()

if (BUILD_FASTSIM)
    message(STATUS "Building FastSim module...")
    add_subdirectory(tang_fastsim)
endif ()

enable_testing()
