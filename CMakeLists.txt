cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(disk_quadrangulation VERSION 0.1 LANGUAGES C CXX)

# Shared libraries
set(BUILD_SHARED_LIBS ON)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
    message(STATUS "CMAKE_BUILD_TYPE set to ${CMAKE_BUILD_TYPE} (default)")
endif()

# Module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/)

set (CMAKE_CXX_STANDARD 17)

# Path to executable directory
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include user options
include(CMakeOptions.txt OPTIONAL)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set (ENABLE_SANITY_CHECKS ON)
  if(ENABLE_SANITY_CHECKS)
    message(STATUS "using sanitizer options")
    add_compile_options("-Wall"
      "-Wextra"
      "-fno-omit-frame-pointer"
      "-Wno-missing-field-initializers"
      "-Wpointer-arith"
      "-Wcast-align"
      "-Wcast-qual"
      "-Wstrict-aliasing"
      "-Wpointer-arith"
      "-Winit-self"
      "-Wredundant-decls"
      "-fsanitize=address"
      "-fsanitize=leak"
      "-fsanitize=undefined")

    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} \
    -fsanitize=address -fsanitize=leak -fsanitize=undefined")
  endif()
endif()

find_package(OpenMP REQUIRED)

# Print useful information
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "Library Type: " ${LIB_TYPE})
MESSAGE(STATUS "Compiler flags:" ${CMAKE_CXX_COMPILE_FLAGS})
MESSAGE(STATUS "Compiler cxx debug flags:" ${CMAKE_CXX_FLAGS_DEBUG})
MESSAGE(STATUS "Compiler cxx release flags:" ${CMAKE_CXX_FLAGS_RELEASE})
MESSAGE(STATUS "Compiler cxx min size flags:" ${CMAKE_CXX_FLAGS_MINSIZEREL})
MESSAGE(STATUS "Compiler cxx flags:" ${CMAKE_CXX_FLAGS})

# fmtlib dependency (formatting and logging)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third_party/fmt)

# gmsh dependency
message(STATUS "----- gmsh library -----")
find_path (GMSH_INCLUDE_DIR
    gmsh.h
    PATHS ${gmsh_API_DIR}
    CMAKE_FIND_ROOT_PATH_BOTH
    )
message(STATUS "gmsh include dir: ${GMSH_INCLUDE_DIR}")
message(STATUS "gmsh lib: ${GMSH_LIBRARY}")
message(STATUS "---------------------------")


# nauty dependency
set(NAUTY_SRC
  third_party/nauty/nauty.c
  third_party/nauty/nautil.c
  third_party/nauty/naugraph.c
  third_party/nauty/schreier.c
  third_party/nauty/naurng.c
  third_party/nauty/nautinv.c
  )
add_library(nauty_dense ${NAUTY_SRC})
target_include_directories(nauty_dense PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/nauty)


set(DQ_SRC_FILES
    src/logging.cpp
    src/geolog.cpp
    src/disk_quadrangulation.cpp
    src/nauty_wrapper.cpp
    )
add_library(dquad_lib ${DQ_SRC_FILES})
target_include_directories(dquad_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(dquad_lib PUBLIC ${GMSH_INCLUDE_DIR})
target_link_libraries(dquad_lib PUBLIC fmt::fmt)
# target_link_libraries(dquad_lib PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(dquad_lib PUBLIC ${GMSH_LIBRARY})
target_link_libraries(dquad_lib PUBLIC nauty_dense)

add_executable(dquad exe/dquad.cpp)
target_link_libraries(dquad dquad_lib)

add_executable(enumerate exe/enumerate.cpp)
target_link_libraries(enumerate dquad_lib)

