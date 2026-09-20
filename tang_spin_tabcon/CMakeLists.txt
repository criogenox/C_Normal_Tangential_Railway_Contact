# Library Option:
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

add_library(kalker_table src/kalker_table.cpp)

target_include_directories(kalker_table PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(tabcon src/main.cpp)

target_include_directories(tabcon PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(tabcon PRIVATE kalker_table)

# Ensure tabcon.dat is available in the build directory for testing
configure_file(tabcon.dat ${CMAKE_CURRENT_BINARY_DIR}/tabcon.dat COPYONLY)


