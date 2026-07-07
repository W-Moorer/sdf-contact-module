# ============================================================================
# SdfOracle Plugin - Standalone and NexDyn subdirectory build
# ============================================================================
# Independent of NexDyn core. Only requires Eigen3.
# Builds:
#   - SdfOracleCore (runtime .sdf IO/query library)
#   - SdfOracle (builder-compatible static library)
#   - SdfOracleDemo (validation executable)
# ============================================================================

cmake_minimum_required(VERSION 3.15)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(SdfOracle VERSION 1.0 LANGUAGES CXX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    add_compile_options(/utf-8 /wd4267 /wd4819 /wd4244 /bigobj)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 CONFIG REQUIRED)
endif()

option(NEXDYN_SDF_ENABLE_OPENMP "Enable OpenMP acceleration in SdfOracle" ON)
if(NEXDYN_SDF_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
endif()

option(NEXDYN_SDF_ENABLE_CUDA "Enable optional CUDA acceleration in SdfOracle" OFF)
if(NEXDYN_SDF_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit QUIET)
    else()
        message(WARNING "NEXDYN_SDF_ENABLE_CUDA=ON but no CUDA compiler was found")
    endif()
endif()

set(SDF_ORACLE_CORE_HEADERS
    include/SdfOracle/SdfCommon.h
    include/SdfOracle/GridSpec.h
    include/SdfOracle/SdfData.h
    include/SdfOracle/SdfOracle.h
    include/SdfOracle/SdfIO.h
    include/SdfOracle/NexDynAdapter.h
)

set(SDF_ORACLE_CORE_SOURCES
    src/SdfCommon.cpp
    src/SdfOracle.cpp
    src/SdfIO.cpp
)

set(SDF_ORACLE_BUILDER_HEADERS
    include/SdfOracle/KdTree.h
    include/SdfOracle/SdfBuilder.h
    include/SdfOracle/SdfSurfaceMesh.h
)

set(SDF_ORACLE_BUILDER_SOURCES
    src/SdfBuilder.cpp
    src/SdfSurfaceMesh.cpp
)

function(configure_sdf_oracle_target target)
    target_include_directories(${target} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(${target} PUBLIC Eigen3::Eigen)
    if(NEXDYN_SDF_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(${target} PUBLIC SDF_ORACLE_HAS_OPENMP=1)
    else()
        target_compile_definitions(${target} PUBLIC SDF_ORACLE_HAS_OPENMP=0)
    endif()
    if(NEXDYN_SDF_ENABLE_CUDA AND CUDAToolkit_FOUND)
        target_link_libraries(${target} PUBLIC CUDA::cudart)
        target_compile_definitions(${target} PUBLIC SDF_ORACLE_ENABLE_CUDA=1)
    else()
        target_compile_definitions(${target} PUBLIC SDF_ORACLE_ENABLE_CUDA=0)
    endif()
endfunction()

add_library(SdfOracleCore STATIC ${SDF_ORACLE_CORE_HEADERS} ${SDF_ORACLE_CORE_SOURCES})
configure_sdf_oracle_target(SdfOracleCore)

add_library(SdfOracle STATIC ${SDF_ORACLE_BUILDER_HEADERS} ${SDF_ORACLE_BUILDER_SOURCES})
configure_sdf_oracle_target(SdfOracle)
target_link_libraries(SdfOracle PUBLIC SdfOracleCore)

add_executable(SdfOracleDemo app/SdfOracleDemo.cpp)
target_link_libraries(SdfOracleDemo PRIVATE SdfOracle)

install(TARGETS SdfOracleCore SdfOracle SdfOracleDemo
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY include/SdfOracle DESTINATION include)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}
    FILES
        ${SDF_ORACLE_CORE_HEADERS}
        ${SDF_ORACLE_CORE_SOURCES}
        ${SDF_ORACLE_BUILDER_HEADERS}
        ${SDF_ORACLE_BUILDER_SOURCES})
source_group("app" FILES app/SdfOracleDemo.cpp)

message(STATUS "SdfOracle plugin configured")
