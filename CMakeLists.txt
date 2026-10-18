cmake_minimum_required(VERSION 3.16)
project(leptonff VERSION 0.27.0 LANGUAGES CXX)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckCXXSourceCompiles)
include(CTest)
find_package(Threads REQUIRED)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(LEPTONFF_BUILD_MATCHED "Build matched W/Z calculations and production records" ON)
option(LEPTONFF_BUILD_TOOLS "Build and install production replay" ON)
option(LEPTONFF_BUILD_REFERENCE_PROGRAMS "Build research/reference programs" OFF)
option(LEPTONFF_BUILD_EXAMPLES "Build introductory scientific examples" ON)
option(LEPTONFF_BUILD_FORTRAN "Build and install the Fortran module (requires a compiler)" OFF)
option(LEPTONFF_BUILD_DOCS "Build the HTML scientific manual with Sphinx" OFF)
if(NOT LEPTONFF_BUILD_MATCHED)
  set(LEPTONFF_BUILD_TOOLS OFF)
  set(LEPTONFF_BUILD_REFERENCE_PROGRAMS OFF)
endif()
foreach(bits IN ITEMS 113 64 53)
  check_cxx_source_compiles("#include <limits>\nstatic_assert(std::numeric_limits<long double>::digits == ${bits}, \"different long double ABI\");\nint main(){}" LEPTONFF_REAL_${bits})
  if(LEPTONFF_REAL_${bits})
    set(LEPTONFF_REAL_MANTISSA ${bits})
    break()
  endif()
endforeach()
if(NOT LEPTONFF_REAL_MANTISSA)
  message(FATAL_ERROR "Unsupported long double ABI")
endif()

# Compile core objects once; the complete archive stays self-contained for
# existing generator relinkers. Consumers choose either core or the full target.
add_library(leptonff_core_objects OBJECT
  src/quadrature.cpp
  src/distribution.cpp
  src/special_functions.cpp
  src/fragmentation.cpp
  src/nlo_lepton.cpp
  src/nlo_qed.cpp
  src/qcd.cpp
  src/hadronic_input.cpp
  src/flavours.cpp
  src/flavour_moment.cpp
  src/flavour_xspace.cpp
  src/moment_evolution.cpp
  src/xspace_evolution.cpp
  src/fft.cpp
  src/evolution.cpp
  src/theory.cpp
  src/fragmentation_c.cpp
  src/api_error.cpp
)
function(leptonff_library_settings target)
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  target_compile_features(${target} PUBLIC cxx_std_17)
  target_compile_definitions(${target} PUBLIC LEPTONFF_REAL_MANTISSA=${LEPTONFF_REAL_MANTISSA})
  set_target_properties(${target} PROPERTIES CXX_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()
leptonff_library_settings(leptonff_core_objects)
target_compile_definitions(leptonff_core_objects PRIVATE LEPTONFF_VERSION="${PROJECT_VERSION}")
add_library(leptonff_core $<TARGET_OBJECTS:leptonff_core_objects>)
leptonff_library_settings(leptonff_core)
set_target_properties(leptonff_core PROPERTIES EXPORT_NAME core VERSION ${PROJECT_VERSION} SOVERSION 0)
add_library(leptonff::core ALIAS leptonff_core)
set(leptonff_install_targets leptonff_core)
if(LEPTONFF_BUILD_MATCHED)
  add_library(leptonff $<TARGET_OBJECTS:leptonff_core_objects>
    src/vector_decay.cpp
    src/hard_dressing.cpp
    src/fiducial_vector.cpp
    src/photon_pair.cpp
    src/pair_moment.cpp
    src/coupled_decay.cpp
    src/mass_spectrum.cpp
    src/lepton_spectrum.cpp
    src/lab_spectrum.cpp
    src/lab_fragmentation.cpp
    src/spin_decay.cpp
    src/photonic_dressing.cpp
    src/photonic_nll.cpp
    src/production_decay.cpp
    src/decay_setup.cpp
    src/response_table.cpp
    src/drell_yan.cpp
    src/c_api.cpp
    src/production_records.cpp
    src/coupled_lab.cpp
    src/vector_current.cpp
    src/decay.cpp
    src/ff_decay.cpp
  )
  leptonff_library_settings(leptonff)
  target_link_libraries(leptonff PUBLIC Threads::Threads)
  # Before 1.0, matched C++ layouts may change between minor releases. Keep
  # older binaries on their original library instead of silently mixing them.
  set_target_properties(leptonff PROPERTIES VERSION ${PROJECT_VERSION}
    SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")
  list(APPEND leptonff_install_targets leptonff)
else()
  add_library(leptonff INTERFACE)
  target_link_libraries(leptonff INTERFACE leptonff::core)
  list(APPEND leptonff_install_targets leptonff)
endif()
add_library(leptonff::leptonff ALIAS leptonff)
install(TARGETS ${leptonff_install_targets} EXPORT leptonffTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/leptonff DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES LICENSE THIRD_PARTY_NOTICES.md DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/leptonff)
install(FILES README.md CITATION.cff CHANGELOG.md DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(DIRECTORY docs/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/docs FILES_MATCHING PATTERN "*.md" PATTERN "*.rst")
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples
  FILES_MATCHING PATTERN "CMakeLists.txt" PATTERN "*.cpp" PATTERN "*.c" PATTERN "*.f90" PATTERN "*.py")
install(DIRECTORY benchmarks/canonical DESTINATION ${CMAKE_INSTALL_DOCDIR}/benchmarks
  FILES_MATCHING PATTERN "*.json" PATTERN "*.md")
install(FILES benchmarks/architecture_0_27.json benchmarks/dyturbo_native_1_4_2.json benchmarks/dyturbo_qed_mass.json
  benchmarks/dyturbo_w_decay.json benchmarks/dyturbo_mass_photos.json
  benchmarks/z_accuracy_scales.json benchmarks/z_radiator_validation.json
  DESTINATION ${CMAKE_INSTALL_DOCDIR}/benchmarks)
install(FILES bindings/fortran/leptonff.f90 DESTINATION ${CMAKE_INSTALL_DATADIR}/leptonff/fortran)
configure_package_config_file(cmake/leptonffConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/leptonffConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/leptonff)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/leptonffConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(EXPORT leptonffTargets NAMESPACE leptonff:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/leptonff)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/leptonffConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/leptonffConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/leptonff)
set(LEPTONFF_PC_LIBRARY leptonff_core)
set(LEPTONFF_PC_THREAD_LIBS "")
if(LEPTONFF_BUILD_MATCHED)
  set(LEPTONFF_PC_LIBRARY leptonff)
  set(LEPTONFF_PC_THREAD_LIBS "${CMAKE_THREAD_LIBS_INIT}")
endif()
file(RELATIVE_PATH LEPTONFF_PC_PREFIX_RELATIVE
  "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig" "${CMAKE_INSTALL_PREFIX}")
configure_file(cmake/leptonff.pc.in leptonff.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/leptonff.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
if(LEPTONFF_BUILD_FORTRAN)
  enable_language(Fortran)
  add_library(leptonff_fortran bindings/fortran/leptonff.f90)
  target_link_libraries(leptonff_fortran PUBLIC leptonff::leptonff)
  set_target_properties(leptonff_fortran PROPERTIES EXPORT_NAME fortran
    Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/fortran-modules)
  target_include_directories(leptonff_fortran INTERFACE
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/fortran-modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/leptonff/fortran/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}>)
  add_library(leptonff::fortran ALIAS leptonff_fortran)
  install(TARGETS leptonff_fortran EXPORT leptonffTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY ${PROJECT_BINARY_DIR}/fortran-modules/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/leptonff/fortran/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION})
endif()
if(LEPTONFF_BUILD_TOOLS)
  add_executable(leptonff_replay apps/replay.cpp)
  target_link_libraries(leptonff_replay PRIVATE leptonff::leptonff)
  if(APPLE)
    set_target_properties(leptonff_replay PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  elseif(UNIX)
    set_target_properties(leptonff_replay PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  endif()
  target_compile_definitions(leptonff_replay PRIVATE LEPTONFF_VERSION="${PROJECT_VERSION}" LEPTONFF_BUILD_TYPE="$<CONFIG>")
  install(TARGETS leptonff_replay RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(PROGRAMS python/replay_production.py DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME leptonff-replay)
endif()
if(LEPTONFF_BUILD_REFERENCE_PROGRAMS)
  add_subdirectory(apps)
endif()
if(LEPTONFF_BUILD_EXAMPLES OR BUILD_TESTING)
  add_subdirectory(examples)
endif()
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
if(LEPTONFF_BUILD_DOCS)
  find_program(SPHINX_BUILD NAMES sphinx-build)
  if(NOT SPHINX_BUILD)
    message(FATAL_ERROR "LEPTONFF_BUILD_DOCS requires sphinx-build; install docs/manual/requirements.txt")
  endif()
  add_custom_target(docs ALL COMMAND ${SPHINX_BUILD} -W --keep-going -b html
    ${PROJECT_SOURCE_DIR}/docs/manual ${PROJECT_BINARY_DIR}/html VERBATIM)
  install(DIRECTORY ${PROJECT_BINARY_DIR}/html/
    DESTINATION ${CMAKE_INSTALL_DOCDIR}/html OPTIONAL)
endif()
