# Reusable helper functions for the Bittree CMake build.
#
# Included once from the root CMakeLists.txt. CMake functions are global, so the
# per-dimension src/ and test/ CMakeLists (added via add_subdirectory) can call
# them. This module is self-contained.

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# --------------------------------------------------------------------------
# bittree_apply_dev_flags(<target>)
#
# Apply the project's strict warning set in Debug builds (GNU/Clang only,
# PRIVATE so it never leaks to consumers), mirroring CXXFLAGS_DEBUG in
# Makefile.site, plus --coverage instrumentation when BITTREE_ENABLE_COVERAGE
# is on. Self-contained: all generator expressions are local.
# (target_link_libraries carries the coverage link flag to stay 3.12-safe;
# target_link_options would require 3.13.)
# --------------------------------------------------------------------------
function(bittree_apply_dev_flags target)
  set(gnuish "$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>")
  set(debug_gnuish "$<AND:$<CONFIG:Debug>,${gnuish}>")

  set(warnings -Wall -Wextra -pedantic -Wconversion)
  if(BITTREE_WERROR)
    list(APPEND warnings -Werror)
  endif()
  target_compile_options(${target} PRIVATE "$<${debug_gnuish}:${warnings}>")

  if(BITTREE_ENABLE_COVERAGE)
    target_compile_options(${target} PRIVATE "$<${gnuish}:--coverage>")
    target_link_libraries(${target}  PRIVATE "$<${gnuish}:--coverage>")
  endif()
endfunction()

# --------------------------------------------------------------------------
# bittree_install_dimension(<target> <header>...)
#
# Install one dimension-specific library as a self-contained CMake package under
# <prefix>/<N>d/ (keyed off BITTREE_CURRENT_DIM), so a downstream
# find_package(Bittree) pointed at <prefix>/<N>d resolves the matching
# Bittree::bittree. The trailing args are the public headers to install
# (relative paths resolve against the calling directory; the generated
# Bittree_constants.h is passed as an absolute path).
# --------------------------------------------------------------------------
function(bittree_install_dimension target)
  set(headers ${ARGN})
  set(dimtag "${BITTREE_CURRENT_DIM}d")

  set(dim_libdir   "${dimtag}/${CMAKE_INSTALL_LIBDIR}")
  set(dim_incdir   "${dimtag}/${CMAKE_INSTALL_INCLUDEDIR}")
  set(dim_bindir   "${dimtag}/${CMAKE_INSTALL_BINDIR}")
  set(dim_cmakedir "${dim_libdir}/cmake/Bittree")

  # The public headers include Bittree_constants.h, so consumers need the
  # installed include dir. _IMPORT_PREFIX resolves to <prefix>, so the
  # interface path is expressed relative to it as <N>d/include.
  target_include_directories(${target} PUBLIC
    $<INSTALL_INTERFACE:${dim_incdir}>)

  install(TARGETS ${target}
          EXPORT bittree_${dimtag}_targets
          ARCHIVE DESTINATION ${dim_libdir}
          LIBRARY DESTINATION ${dim_libdir}
          RUNTIME DESTINATION ${dim_bindir}
          INCLUDES DESTINATION ${dim_incdir})

  install(FILES ${headers} DESTINATION ${dim_incdir})

  install(EXPORT bittree_${dimtag}_targets
          NAMESPACE Bittree::
          FILE BittreeTargets.cmake
          DESTINATION ${dim_cmakedir})

  # Package config + version file, generated into a per-dim build subdir (the
  # subdir binary dir is already unique per dimension) so they don't collide.
  set(cfgbuild "${CMAKE_CURRENT_BINARY_DIR}/cmake-pkg")
  set(BITTREE_DIM ${BITTREE_CURRENT_DIM})   # consumed by BittreeConfig.cmake.in
  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/BittreeConfig.cmake.in"
    "${cfgbuild}/BittreeConfig.cmake"
    INSTALL_DESTINATION ${dim_cmakedir})
  write_basic_package_version_file(
    "${cfgbuild}/BittreeConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
  install(FILES
          "${cfgbuild}/BittreeConfig.cmake"
          "${cfgbuild}/BittreeConfigVersion.cmake"
          DESTINATION ${dim_cmakedir})
endfunction()
