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

# Bittree library for one dimension.
#
# Driven by the root CMakeLists.txt, which sets BITTREE_CURRENT_DIM and calls
# add_subdirectory(src ...) once per dimension with a distinct binary dir.
# This file is not meant to be added standalone.

set(dimtag "${BITTREE_CURRENT_DIM}d")
set(libtgt "bittree_${dimtag}")

# Generated dimension header (equivalent to what setup.py writes); consumed by
# cmake/Bittree_constants.h.in via @BITTREE_DIM@. The binary dir is already
# unique per dimension (build/src/<N>d), so this never collides.
set(gendir "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(BITTREE_DIM ${BITTREE_CURRENT_DIM})
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/Bittree_constants.h.in"
  "${gendir}/Bittree_constants.h"
  @ONLY)

add_library(${libtgt} STATIC
  Bittree_BitArray.cpp
  Bittree_MortonTree.cpp
  Bittree_BittreeAmr.cpp
  Bittree_fi.cpp)

set_target_properties(${libtgt} PROPERTIES
  OUTPUT_NAME bittree     # -> libbittree.a
  EXPORT_NAME bittree)    # -> Bittree::bittree

target_compile_features(${libtgt} PUBLIC cxx_std_17)

# PUBLIC: a public header (Bittree_MortonTree.h) includes Bittree_constants.h,
# so consumers need both the src headers and the generated one. The matching
# INSTALL_INTERFACE include is added by bittree_install_dimension().
target_include_directories(${libtgt} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${gendir}>)

target_link_libraries(${libtgt} PUBLIC MPI::MPI_CXX)

if(BITTREE_ALT_MORTON_ORDER)
  target_compile_definitions(${libtgt} PUBLIC ALT_MORTON_ORDER)
endif()

bittree_apply_dev_flags(${libtgt})

# Install this dimension as a self-contained package under <prefix>/<N>d/.
bittree_install_dimension(${libtgt}
  Bittree_BitArray.h
  Bittree_Bits.h
  Bittree_BittreeAmr.h
  Bittree_MortonTree.h
  Bittree_Prelude.h
  Bittree_fi.h
  "${gendir}/Bittree_constants.h")
