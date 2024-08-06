cmake_minimum_required (VERSION 3.1)
project(SBWT)
include(ExternalProject)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  option(GIT_SUBMODULE "Check submodules during build" ON)
  if(GIT_SUBMODULE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE SBWT_BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message(STATUS "SBWT version: ${SBWT_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(SBWT_BUILD_VERSION "unknown")
endif()


## Set maximum k-mer length
if(NOT MAX_KMER_LENGTH)
  set(MAX_KMER_LENGTH 255)
endif()
if(MAX_KMER_LENGTH GREATER 255)
  message(${MAX_KMER_LENGTH})
  message(FATAL_ERROR "Error: MAX_KMER_LENGTH can be at most 255")
endif()
add_definitions(-DMAX_KMER_LENGTH=${MAX_KMER_LENGTH}) # Define for compiler.
message(STATUS "Maximum k-mer length is: " ${MAX_KMER_LENGTH})
message(STATUS "To change the maximum k-mer length, use the option -DMAX_KMER_LENGTH=n")
message(STATUS "(n can be up to 255)")

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
  # set(CMAKE_BUILD_TYPE Debug)
endif()

# if(CMAKE_BUILD_TYPE MATCHES Release)
#   add_compile_options(/O1)
# endif()

# if(CMAKE_BUILD_TYPE MATCHES Debug)
#   add_compile_options(/W4 )
# endif()

if(CMAKE_BUILD_TYPE MATCHES Release)
  add_compile_options(-O3)
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g -Wall -Wextra -Wno-sign-compare)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

## Require c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fix for GCC (see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90050 )
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")
endif()

# if(MSVC)
#   set(CMAKE_CXX_FLAGS_RELEASE "/MT")
#   set(CMAKE_CXX_FLAGS_DEBUG "/MTd")
# endif()

# find_package(Threads REQUIRED)

## Add local dependencies as targets
add_subdirectory(${PROJECT_SOURCE_DIR}/sdsl-lite
  ${PROJECT_BINARY_DIR}/external/sdsl-lite/build)

# add_subdirectory(${PROJECT_SOURCE_DIR}/KMC
#   ${PROJECT_BINARY_DIR}/external/KMC/build)

# Our sources, excluding the file with the main function
# set(SBWT_alphabet_SOURCES
#   sbwt_alphabet.cpp
# )

add_executable(sbwt
  main.cpp
  ${SBWT_alphabet_SOURCES})

  target_include_directories(sbwt PRIVATE
  ${PROJECT_SOURCE_DIR}/include
#   ${PROJECT_SOURCE_DIR}/include/sbwt/libwheeler
  ${PROJECT_SOURCE_DIR}/sdsl-lite/include
#   ${PROJECT_BINARY_DIR}/external/sdsl-lite/build/external/libdivsufsort/include
  ${PROJECT_BINARY_DIR}/include
  ${PROJECT_SOURCE_DIR}
)
target_link_libraries(sbwt PRIVATE
  sdsl
#   ${BZIP2}
#   ${ZLIB}
#   Threads::Threads
  ${CXX_FILESYSTEM_LIBRARIES}
)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sbwt PUBLIC OpenMP::OpenMP_CXX)
endif()



# add_library(sbwt_static STATIC ${SBWT_SOURCES})
# target_include_directories(sbwt_static PRIVATE
#   ${PROJECT_SOURCE_DIR}/include/sbwt/
#   ${PROJECT_SOURCE_DIR}/include/sbwt/libwheeler
#   ${PROJECT_SOURCE_DIR}/sdsl-lite/include
#   ${PROJECT_BINARY_DIR}/external/sdsl-lite/build/external/libdivsufsort/include
#   ${PROJECT_BINARY_DIR}/include
#   ${PROJECT_SOURCE_DIR}
#   ${PROJECT_SOURCE_DIR}/KMC/include)
# target_link_libraries(sbwt_static PRIVATE
#   sdsl
#   ${BZIP2}
#   ${ZLIB}
#   Threads::Threads
#   ${CXX_FILESYSTEM_LIBRARIES}
#   kmc_tools
#   kmc_core)
