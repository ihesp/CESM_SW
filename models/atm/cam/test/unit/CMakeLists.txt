cmake_minimum_required(VERSION 2.8)
project(CAM_tests Fortran C)

# Unit test boilerplate.
list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)

# Useful paths.
set(CAMROOT ../../)
set(CESMROOT ${CAMROOT}../../../)
set(SHAREROOT ${CESMROOT}/models/csm_share)

# Add share code to source list.
include_directories(${SHAREROOT}/include ${SHAREROOT}/shr)
add_subdirectory(${SHAREROOT}/shr csm_share)

# CAM physics.
add_subdirectory(${CAMROOT}src/physics/cam physics_cam)
add_subdirectory(${CAMROOT}src/utils utils_cam)

# Stub/mock modules.
add_subdirectory(mock)

# Build part of csm_share as a library.
add_library(csm_share ${share_sources})
get_target_property(includes csm_share INCLUDE_DIRECTORIES)
if(NOT includes)
  unset(includes)
endif()
list(APPEND includes "${PFUNIT_INCLUDE_DIRS}")
set_target_properties(csm_share PROPERTIES
  INCLUDE_DIRECTORIES "${includes}")

# Pick up csm_share module files.
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Test subdirectories.
add_subdirectory(coords_1d)
add_subdirectory(linear_1d_operators)
add_subdirectory(micro_mg_data)
add_subdirectory(vdiff_lu_solver)
