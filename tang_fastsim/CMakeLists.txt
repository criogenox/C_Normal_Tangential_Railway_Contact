# ================== BUILD OPTIONS (switch) ==================
option(BUILD_STATIC_LIB "Build the static library FastSim" ON)
option(BUILD_TESTS "Build GoogleTest suite" ON)

# Paths
set(LIBS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libs")
set(HEADERS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(UTILS_PATH "${CMAKE_CURRENT_SOURCE_DIR}/utils")
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(TEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tests")

include_directories(${HEADERS_PATH} ${UTILS_PATH})

# ====================== STATIC LIBRARY ======================
if (BUILD_STATIC_LIB)
    add_library(FastSim_lib STATIC ${SRC_PATH}/fastsim.cpp)

    set_target_properties(FastSim_lib PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY ${LIBS_PATH}
            OUTPUT_NAME "fastsim"
            EXPORT_NAME FastSim
    )

    # Install library
    install(TARGETS FastSim_lib
            EXPORT FastSimTargets
            ARCHIVE DESTINATION ${LIBS_PATH}
    )
endif ()

# ====================== GOOGLE TESTS ======================
if (BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(googletest
            GIT_REPOSITORY https://github.com/google/googletest
            GIT_TAG v1.14.0)

    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

    FetchContent_MakeAvailable(googletest)

    file(GLOB_RECURSE TEST_SOURCES "${TEST_PATH}/*.cpp")
    add_executable(FastSimTest ${TEST_SOURCES})

    target_link_libraries(FastSimTest PRIVATE FastSim_lib gtest_main)

    include(GoogleTest)
    enable_testing()
    gtest_discover_tests(FastSimTest)
endif ()