cmake_minimum_required(VERSION 3.16)

# ScalFMM 1.5 lite: sequential/OpenMP FMM on a fixed set of particles,
# Rotation (spherical harmonics) and Chebyshev (interpolation) kernels only.
project(SCALFMM VERSION 1.5.2 LANGUAGES CXX)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR "ScalFMM requires an out of source build (mkdir Build && cd Build && cmake ..)")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Release, Debug, RelWithDebInfo)" FORCE)
endif()

#===========================================================================
# Options
#===========================================================================
option(SCALFMM_USE_ASSERT    "Enable safety tests during execution"                ON)
option(SCALFMM_USE_LOG       "Enable the debug log (FLog)"                         OFF)
option(SCALFMM_USE_MEM_STATS "Trace memory allocations"                            OFF)
option(SCALFMM_USE_SSE       "Use SSE intrinsics in P2P"                           OFF)
option(SCALFMM_USE_AVX       "Use AVX intrinsics in P2P"                           ON)
option(SCALFMM_USE_NATIVE    "Compile with -march=native"                          ON)
option(SCALFMM_BUILD_TESTS   "Build the Rotation/Chebyshev validation tests"       ON)
set(SCALFMM_BLAS_MANGLING "ADD_" CACHE STRING "Fortran mangling of BLAS/LAPACK symbols (ADD_, UPCASE, NOCHANGE)")
set_property(CACHE SCALFMM_BLAS_MANGLING PROPERTY STRINGS ADD_ UPCASE NOCHANGE)

if(SCALFMM_USE_SSE AND SCALFMM_USE_AVX)
  message(FATAL_ERROR "SCALFMM_USE_SSE and SCALFMM_USE_AVX are mutually exclusive")
endif()

#===========================================================================
# Dependencies: OpenMP, BLAS + LAPACK (Chebyshev kernel: SVD, gemm...)
#===========================================================================
find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(SCALFMM_USE_BLAS ON)
set(SCALFMM_BLAS_ADD_ OFF)
set(SCALFMM_BLAS_UPCASE OFF)
set(SCALFMM_BLAS_NOCHANGE OFF)
set(SCALFMM_BLAS_${SCALFMM_BLAS_MANGLING} ON)

#===========================================================================
# Compile flags
#===========================================================================
set(SCALFMM_CXX_FLAGS -m64 -fpic -ffast-math -funroll-loops -ftree-vectorize)
if(SCALFMM_USE_NATIVE)
  list(APPEND SCALFMM_CXX_FLAGS -march=native -mtune=native)
elseif(SCALFMM_USE_AVX)
  list(APPEND SCALFMM_CXX_FLAGS -mavx)
elseif(SCALFMM_USE_SSE)
  list(APPEND SCALFMM_CXX_FLAGS -msse4)
endif()

string(REPLACE ";" " " SCALFMM_COMPILE_FLAGS "${SCALFMM_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
string(REPLACE ";" " " SCALFMM_COMPILE_LIBS  "${LAPACK_LIBRARIES};${BLAS_LIBRARIES}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/Src/ScalFmmConfig.h)

#===========================================================================
# Library
#===========================================================================
add_library(scalfmm STATIC
  Src/Utils/FAlgorithmTimers.cpp
  Src/Utils/FLog.cpp
  Src/Utils/FMemStats.cpp)
set_target_properties(scalfmm PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/${CMAKE_BUILD_TYPE})
target_compile_features(scalfmm PUBLIC cxx_std_14)
target_compile_options(scalfmm PUBLIC ${SCALFMM_CXX_FLAGS})
target_include_directories(scalfmm PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/Src>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(scalfmm PUBLIC OpenMP::OpenMP_CXX ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

#===========================================================================
# Tests
#===========================================================================
if(SCALFMM_BUILD_TESTS)
  enable_testing()
  foreach(t testRotationTsm testChebyshevTsm)
    add_executable(${t} Tests/${t}.cpp)
    target_link_libraries(${t} PRIVATE scalfmm)
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
endif()

#===========================================================================
# Install (same layout as ScalFMM 1.5: <prefix>/include/<dir>, <prefix>/lib)
#===========================================================================
install(TARGETS scalfmm EXPORT ScalFMMTargets ARCHIVE DESTINATION lib)
install(DIRECTORY Src/ DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Src/ScalFmmConfig.h DESTINATION include)
install(EXPORT ScalFMMTargets NAMESPACE scalfmm:: DESTINATION lib/cmake/ScalFMM)
