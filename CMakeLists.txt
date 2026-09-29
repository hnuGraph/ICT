cmake_minimum_required(VERSION 3.10)
project(DSQL_META)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)

add_library(project_options INTERFACE)

if(MSVC)
    target_compile_options(project_options INTERFACE /O2 /W4)
else()
    target_compile_options(project_options INTERFACE
        "$<$<CONFIG:Release>:-O3;-march=native>"
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

include_directories(${CMAKE_SOURCE_DIR})

# Utility
set(UTIL_SRC
    utils/globals.cpp
)
add_library(utility STATIC ${UTIL_SRC})
target_link_libraries(utility PUBLIC project_options)

# Graph
set(GRAPH_SRC
    graph/graph.cpp
)
add_library(graph STATIC ${GRAPH_SRC})
target_link_libraries(graph PUBLIC project_options)

# Matching
set(MATCHING_SRC
    matching/matching.cpp
    matching/DsqlMeta.cpp
)
add_library(matching STATIC ${MATCHING_SRC})
target_link_libraries(matching PUBLIC project_options)

# Main executable
add_executable(DSQL_META main.cpp)

target_link_libraries(
    DSQL_META
    PRIVATE
        utility
        graph
        matching
        Threads::Threads
)