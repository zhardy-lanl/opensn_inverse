cmake_minimum_required(VERSION 3.12)

project(inverse_problem
    VERSION 0.0.1
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CheckTypeSize)
include(CheckSymbolExists)

# dependencies
find_package(MPI REQUIRED)
find_package(OpenSn REQUIRED)
find_package(Lua 5.4 REQUIRED)

# libinverse
file(GLOB_RECURSE LIBINVERSE_SRCS CONFIGURE_DEPENDS inverse_solver/*.cc)
add_library(libinverse OBJECT ${LIBINVERSE_SRCS})

target_include_directories(libinverse
    PRIVATE
    $<INSTALL_INTERFACE:include/inverse>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}
    ${OPENSN_INCLUDE_DIR}
)

target_link_libraries(libinverse
    PUBLIC
    opensn::libopensn
    opensn::libopensnlua
    PRIVATE
    MPI::MPI_CXX
)

# opensn_inverse binary
add_executable(inverse "main.cc")

target_include_directories(inverse
    PRIVATE
    $<INSTALL_INTERFACE:include/inverse>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    ${OPENSN_INCLUDE_DIR}
)

target_link_libraries(inverse
    PUBLIC
    opensn::libopensn
    opensn::libopensnlua
    libinverse
    PRIVATE
    MPI::MPI_CXX
)
