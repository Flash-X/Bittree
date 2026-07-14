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

# Bittree GoogleTest unit test for one dimension.
#
# Driven by the root CMakeLists.txt (BITTREE_CURRENT_DIM + add_subdirectory per
# dimension). The library target bittree_<N>d and BITTREE_GTEST_LIBS come from
# the root/src scope.

set(dimtag "${BITTREE_CURRENT_DIM}d")
set(testtgt "bittree_unit_test_${dimtag}")

add_executable(${testtgt}
  test.cpp
  main.cpp)

set_target_properties(${testtgt} PROPERTIES OUTPUT_NAME bittree_unit_test.x)

target_include_directories(${testtgt} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})   # macros.h
target_link_libraries(${testtgt} PRIVATE bittree_${dimtag} ${BITTREE_GTEST_LIBS})

bittree_apply_dev_flags(${testtgt})

# Run the single binary directly (it does its own singleton MPI_Init), mirroring
# the Make build's `make test`. Deliberately not gtest_discover_tests, which
# would run the MPI binary at build time (bad on HPC login nodes).
add_test(NAME ${testtgt} COMMAND $<TARGET_FILE:${testtgt}>)
