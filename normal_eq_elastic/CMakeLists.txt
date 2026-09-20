find_package(OpenMP REQUIRED)

add_executable(simrail_contact
        src/main.cpp
        src/utils/MathUtils.cpp
        src/contact/EquivalentContact.cpp
        src/geometry/Interpolator.cpp
        src/geometry/ProfileLoader.cpp
)

# Copy profile files to build directory
add_custom_command(TARGET simrail_contact POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/S1002.txt
        ${CMAKE_CURRENT_SOURCE_DIR}/UIC60.txt
        $<TARGET_FILE_DIR:simrail_contact>
)

target_include_directories(simrail_contact PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/contact
        ${CMAKE_CURRENT_SOURCE_DIR}/include/geometry
        ${CMAKE_CURRENT_SOURCE_DIR}/include/utils
)

target_link_libraries(simrail_contact PRIVATE OpenMP::OpenMP_CXX)

# Optional: Enable OpenMP in compile definitions
target_compile_definitions(simrail_contact PRIVATE USE_OPENMP)