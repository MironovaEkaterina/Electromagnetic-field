cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME Electromagnetic_field)
project(${PROJECT_NAME})

# TODO(Kornyakov): not sure if these lines are needed
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(INSTALL_DIR "${CMAKE_CURRENT_LIST_DIR}/3rdparty")
file(MAKE_DIRECTORY "${INSTALL_DIR}")

set(FFTW_VERSION 3.3.8)
message(STATUS "Downloading fftw-${FFTW_VERSION}")
include(ExternalProject)
ExternalProject_Add(project_fftw
URL "http://fftw.org/fftw-${FFTW_VERSION}.tar.gz"
PREFIX ${INSTALL_DIR}/fftw/src
CMAKE_ARGS
"-DCMAKE_CFLAGS=${CMAKE_C_FLAGS}"
"-DBUILD_SHARED_LIBS=OFF"
"-DENABLE_AVX2=ON"
"-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}/fftw"
"-DENABLE_OPENMP=ON"
"-DBUILD_TESTS=OFF"
)
install(DIRECTORY "${INSTALL_DIR}" DESTINATION .)

set(FFTW3_LIB ${INSTALL_DIR}/fftw/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX})
set(FFTW3_OMP_LIB ${INSTALL_DIR}/fftw/lib/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3_omp${CMAKE_STATIC_LIBRARY_SUFFIX})

set(FFT_INCLUDES ${INSTALL_DIR}/fftw/include)
include_directories("${FFT_INCLUDES}")
set(FFT_LIBRARIES ${FFTW3_OMP_LIB} ${FFTW3_LIB})

add_library(fftw3 STATIC IMPORTED)
set_target_properties(fftw3 PROPERTIES IMPORTED_LOCATION "${FFTW3_LIB}")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})

#set(MP2_LIBRARY "${PROJECT_NAME}")
set(MP2_CUSTOM "${PROJECT_NAME}")
set(MP2_TESTS   "test_${PROJECT_NAME}")
set(MP2_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")

include_directories("${MP2_INCLUDE}" gtest)

# BUILD
add_subdirectory(include)
#add_subdirectory(src)
add_subdirectory(samples)
add_subdirectory(gtest)
add_subdirectory(test)

# REPORT
message( STATUS "")
message( STATUS "General configuration for ${PROJECT_NAME}")
message( STATUS "======================================")
message( STATUS "")
message( STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message( STATUS "")

