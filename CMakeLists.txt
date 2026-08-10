# Copyright 2026 UChicago Argonne, LLC and contributors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.12...3.28)

project(Bittree
        VERSION 1.0.0
        DESCRIPTION "Memory-efficient bittree data structure for AMR octree meshes"
        LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --------------------------------------------------------------------------
# Project-wide defaults
#
# Bittree may be consumed as a subproject (add_subdirectory / FetchContent).
# Anything with a global effect (cache writes, enable_testing) is therefore
# applied only when Bittree is the top-level project, so it never hijacks a
# parent project's settings. PROJECT_IS_TOP_LEVEL is CMake 3.21+, so derive it.
# --------------------------------------------------------------------------
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(BITTREE_IS_TOP_LEVEL TRUE)
else()
  set(BITTREE_IS_TOP_LEVEL FALSE)
endif()

# Directory-scoped: affects only Bittree's own targets, safe as a subproject.
set(CMAKE_CXX_EXTENSIONS OFF)   # emit -std=c++17, not gnu++17 (matches the Make build)

if(BITTREE_IS_TOP_LEVEL)
  if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed; use e.g. `cmake -S . -B build`.")
  endif()

  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  # Default to an optimized build (Release ~ prod); pass -DCMAKE_BUILD_TYPE=Debug
  # to get the strict warnings, mirroring `setup.py ... --debug`.
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
  endif()

  # Default the install prefix inside the build tree (unless given), so installs
  # are self-contained and never need a system prefix / elevated privileges.
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install prefix" FORCE)
    message(STATUS "CMAKE_INSTALL_PREFIX defaulted to '${CMAKE_INSTALL_PREFIX}'.")
  endif()
endif()

# --------------------------------------------------------------------------
# Options (mirror the setup.py switches)
# --------------------------------------------------------------------------
set(BITTREE_NDIM "1;2;3" CACHE STRING
    "Semicolon-separated list of dimensions to build/install; each must be 1, 2, or 3.")
# Tests are on by default for a standalone build, off when embedded in a parent
# project (a consumer should not inherit Bittree's test suite or its gtest dep).
option(BITTREE_BUILD_TESTS      "Build the GoogleTest unit test (setup.py test mode)." ${BITTREE_IS_TOP_LEVEL})
option(BITTREE_ENABLE_COVERAGE  "Enable code-coverage instrumentation (setup.py --coverage)." OFF)
option(BITTREE_WERROR           "Treat warnings as errors in Debug builds." ON)
option(BITTREE_ALT_MORTON_ORDER "Define ALT_MORTON_ORDER." OFF)

list(REMOVE_DUPLICATES BITTREE_NDIM)
foreach(_dim IN LISTS BITTREE_NDIM)
  if(NOT _dim MATCHES "^[123]$")
    message(FATAL_ERROR
            "BITTREE_NDIM entries must each be 1, 2, or 3 (got '${_dim}' in '${BITTREE_NDIM}').")
  endif()
endforeach()
message(STATUS "Bittree: building dimension(s) [${BITTREE_NDIM}]")

include(BittreeHelpers)

# --------------------------------------------------------------------------
# Dependencies (found once, shared by every dimension)
# --------------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS CXX)

if(BITTREE_BUILD_TESTS)
  if(BITTREE_IS_TOP_LEVEL)
    enable_testing()
  endif()
  find_package(GTest REQUIRED)
  # Homebrew/upstream export GTest::gtest_main; CMake's FindGTest module
  # exports GTest::Main. The "*_main" target pulls in the core gtest library.
  if(TARGET GTest::gtest_main)
    set(BITTREE_GTEST_LIBS GTest::gtest_main)
  else()
    set(BITTREE_GTEST_LIBS GTest::Main)
  endif()
endif()

# --------------------------------------------------------------------------
# Build each requested dimension. add_subdirectory is called once per dimension
# with a distinct binary dir; the current dimension is passed down via
# BITTREE_CURRENT_DIM (consumed by src/ and test/ CMakeLists).
# --------------------------------------------------------------------------
# CMAKE_CURRENT_BINARY_DIR (not CMAKE_BINARY_DIR): as a subproject the latter is
# the *parent's* build root, which would scatter Bittree's build tree into it.
foreach(_dim IN LISTS BITTREE_NDIM)
  set(BITTREE_CURRENT_DIM ${_dim})
  add_subdirectory(src "${CMAKE_CURRENT_BINARY_DIR}/src/${_dim}d")
  if(BITTREE_BUILD_TESTS)
    add_subdirectory(test "${CMAKE_CURRENT_BINARY_DIR}/test/${_dim}d")
  endif()
endforeach()

# --------------------------------------------------------------------------
# Coverage report target (mirrors the Makefile `coverage` rule); aggregates
# every configured dimension in the build tree.
# --------------------------------------------------------------------------
if(BITTREE_ENABLE_COVERAGE)
  find_program(LCOV_PROGRAM lcov)
  find_program(GENHTML_PROGRAM genhtml)
  if(LCOV_PROGRAM AND GENHTML_PROGRAM)
    add_custom_target(coverage
      COMMAND ${LCOV_PROGRAM} -o lcov_temp.info -c -d . -b "${CMAKE_CURRENT_SOURCE_DIR}" --no-external
      COMMAND ${GENHTML_PROGRAM} -o Coverage_Report lcov_temp.info
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating coverage report with lcov/genhtml")
  else()
    message(WARNING "BITTREE_ENABLE_COVERAGE is ON but lcov/genhtml were not found; "
      "the 'coverage' target will be unavailable.")
  endif()
endif()
