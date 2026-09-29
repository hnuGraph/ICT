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
        -Wall -Wextra -Wpedantic
    )
endif()

include_directories(${CMAKE_SOURCE_DIR})

#utility.o
set(UTIL_SRC 
    utils/globals.cpp
    )
add_library(utility STATIC ${UTIL_SRC})
target_link_libraries(utility PUBLIC project_options)

#graph.o
set(GRAPH_SRC
    graph/graph.cpp
    )
add_library(graph STATIC ${GRAPH_SRC})
target_link_libraries(graph PUBLIC project_options)


#matching.o
set(MATCHING_SRC
    matching/matching.cpp
    matching/DsqlMeta.cpp
    )
add_library(matching STATIC ${MATCHING_SRC})
target_link_libraries(matching PUBLIC project_options)


add_executable(DSQL_META main.cpp)
target_link_libraries(DSQL_META PRIVATE utility graph matching Threads::Threads)

include(CTest)
if(BUILD_TESTING)
    add_executable(ict_dynamic_coverage_regression ../../tests/dynamic_coverage_regression.cpp)
    target_compile_definitions(ict_dynamic_coverage_regression PRIVATE USE_ICT)
    target_link_libraries(ict_dynamic_coverage_regression PRIVATE utility graph matching Threads::Threads project_options)
    add_test(NAME ict_dynamic_coverage_regression COMMAND ict_dynamic_coverage_regression
        ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/dynamic_coverage/data.graph
        ${CMAKE_CURRENT_SOURCE_DIR}/../../tests/dynamic_coverage/query.graph)
    add_executable(ict_cache_test tests/ict_cache_test.cpp)
    target_link_libraries(ict_cache_test PRIVATE graph project_options)
    add_test(
        NAME ict_cache_test
        COMMAND ict_cache_test ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/data.graph
    )

    add_executable(ict_coverage_smoke tests/coverage_driver.cpp)
    target_compile_definitions(ict_coverage_smoke PRIVATE USE_ICT)
    target_link_libraries(
        ict_coverage_smoke PRIVATE utility graph matching Threads::Threads project_options
    )
    add_test(
        NAME ict_coverage_smoke_q1
        COMMAND ict_coverage_smoke
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/data.graph
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/smoke/queries/q1.graph
            3
    )

endif()
