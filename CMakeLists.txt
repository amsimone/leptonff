cmake_minimum_required(VERSION 3.16)
project(leptonff VERSION 0.25.0 LANGUAGES CXX)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckCXXSourceCompiles)
option(BUILD_SHARED_LIBS "Build the leptonff shared library" OFF)
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

add_library(leptonff
  src/quadrature.cpp
  src/distribution.cpp
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
  src/drell_yan.cpp
  src/c_api.cpp
  src/production_records.cpp
  src/coupled_lab.cpp
  src/vector_current.cpp
  src/decay.cpp)
add_library(leptonff::leptonff ALIAS leptonff)
target_include_directories(leptonff PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(leptonff PUBLIC cxx_std_17)
target_compile_definitions(leptonff PUBLIC LEPTONFF_REAL_MANTISSA=${LEPTONFF_REAL_MANTISSA})
set_target_properties(leptonff PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(leptonff PROPERTIES POSITION_INDEPENDENT_CODE ON
  VERSION ${PROJECT_VERSION} SOVERSION 0)
install(TARGETS leptonff EXPORT leptonffTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/leptonff DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES LICENSE THIRD_PARTY_NOTICES.md
  DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/leptonff)
install(FILES README.md CITATION.cff CHANGELOG.md
  DESTINATION ${CMAKE_INSTALL_DOCDIR})
install(DIRECTORY docs/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/docs
  FILES_MATCHING PATTERN "*.md")
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DOCDIR}/examples
  FILES_MATCHING PATTERN "CMakeLists.txt" PATTERN "*.cpp" PATTERN "*.c" PATTERN "*.f90")
configure_package_config_file(cmake/leptonffConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/leptonffConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/leptonff)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/leptonffConfigVersion.cmake
  VERSION ${PROJECT_VERSION} COMPATIBILITY SameMinorVersion)
install(EXPORT leptonffTargets NAMESPACE leptonff:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/leptonff)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/leptonffConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/leptonffConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/leptonff)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(leptonff PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(leptonff_replay apps/replay.cpp)
if(APPLE)
  set_target_properties(leptonff_replay PROPERTIES INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
  set_target_properties(leptonff_replay PROPERTIES INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
install(TARGETS leptonff_replay RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(PROGRAMS python/replay_production.py DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME leptonff-replay)
add_executable(leptonff_photonic_nll apps/photonic_nll_reference.cpp)
add_executable(leptonff_hadronic apps/hadronic_reference.cpp)
add_executable(leptonff_flavours apps/flavour_reference.cpp)
add_executable(leptonff_spin apps/spin_decay_reference.cpp)
add_executable(leptonff_decay apps/decay_reference.cpp)
add_executable(leptonff_ll apps/ll_reference.cpp)
add_executable(leptonff_match apps/matching_reference.cpp)
add_executable(leptonff_hard apps/vector_mass_reference.cpp)
add_executable(leptonff_fiducial apps/fiducial_reference.cpp)
add_executable(leptonff_dressed apps/hard_dressed_reference.cpp)
add_executable(leptonff_kernels apps/nlo_kernel_reference.cpp)
add_executable(leptonff_moments apps/moment_reference.cpp)
add_executable(leptonff_xspace apps/xspace_reference.cpp)
add_executable(leptonff_photon_pair apps/photon_pair_reference.cpp)
add_executable(leptonff_pair_moments apps/pair_moment_reference.cpp)
add_executable(leptonff_coupled apps/coupled_decay_reference.cpp)
add_executable(leptonff_mass_spectrum apps/mass_spectrum_reference.cpp)
add_executable(leptonff_lepton_spectrum apps/lepton_spectrum_reference.cpp)
foreach(app IN ITEMS leptonff_replay leptonff_photonic_nll leptonff_hadronic leptonff_flavours leptonff_spin leptonff_lepton_spectrum leptonff_mass_spectrum leptonff_coupled leptonff_pair_moments leptonff_photon_pair leptonff_xspace leptonff_moments leptonff_kernels leptonff_dressed leptonff_decay leptonff_ll leptonff_match leptonff_hard leptonff_fiducial)
  target_link_libraries(${app} PRIVATE leptonff)
  target_compile_definitions(${app} PRIVATE
    LEPTONFF_VERSION="${PROJECT_VERSION}" LEPTONFF_BUILD_TYPE="$<CONFIG>")
endforeach()

include(CTest)
if(BUILD_TESTING)
  enable_language(C)
  add_executable(test_c_api tests/test_c_api.c)
  target_link_libraries(test_c_api PRIVATE leptonff)
  if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_c_api PRIVATE -UNDEBUG)
  endif()
  add_test(NAME c_api COMMAND test_c_api)
  foreach(name IN ITEMS drell_yan production_records production_decay photonic_nll qcd_hadronic flavours flavour_evolution lab_spectrum fraction_cells coupled_lab photonic_dressing spin_decay lab_fragmentation endpoint_geometry lepton_spectrum mass_spectrum coupled_decay pair_moment photon_pair distributions fragmentation nlo_lepton nlo_qed moment_evolution xspace_evolution decay evolution matching vector_decay fiducial_vector hard_dressing)
    add_executable(test_${name} tests/test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE leptonff)
    add_test(NAME ${name} COMMAND test_${name})
  endforeach()
  add_test(NAME spin_decay_reference COMMAND leptonff_spin)
  add_test(NAME nlo_kernel_reference COMMAND leptonff_kernels)
  add_test(NAME photon_pair_reference COMMAND leptonff_photon_pair
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  add_test(NAME pair_moment_reference COMMAND leptonff_pair_moments
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  add_test(NAME coupled_decay_reference COMMAND leptonff_coupled
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf 256 2 --dump)
  add_test(NAME lepton_spectrum_reference COMMAND leptonff_lepton_spectrum
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf ${CMAKE_CURRENT_SOURCE_DIR}/configs/lepton_bins.conf 256 1 --dump)
  add_test(NAME mass_spectrum_reference COMMAND leptonff_mass_spectrum
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf ${CMAKE_CURRENT_SOURCE_DIR}/configs/mass_bins.conf 256 1 --dump)
  add_test(NAME moment_reference COMMAND leptonff_moments)
  add_test(NAME xspace_reference COMMAND leptonff_xspace
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf 64)
  add_test(NAME decay_reference COMMAND leptonff_decay
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  add_test(NAME ll_reference COMMAND leptonff_ll
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  add_test(NAME matching_reference COMMAND leptonff_match
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  add_test(NAME vector_mass_reference COMMAND leptonff_hard
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/vector_mass.conf)
  add_test(NAME vector_mass_reject_angular_cuts COMMAND leptonff_hard
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  set_tests_properties(vector_mass_reject_angular_cuts PROPERTIES WILL_FAIL TRUE)
  add_test(NAME fiducial_reference COMMAND leptonff_fiducial
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
  add_test(NAME hard_dressed_reference COMMAND leptonff_dressed
    ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf 128 0.15 32 4)
  find_package(Python3 3.10 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME production_replay_cli COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_replay_cli.py $<TARGET_FILE:leptonff_replay>)
    add_test(NAME dyturbo_source_preparation COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_dyturbo_prepare.py)
    find_program(LEPTONFF_FORTRAN_EXAMPLE_COMPILER NAMES gfortran-16 gfortran)
    set(LEPTONFF_INSTALL_TEST_ARGUMENTS "")
    if(LEPTONFF_FORTRAN_EXAMPLE_COMPILER)
      list(APPEND LEPTONFF_INSTALL_TEST_ARGUMENTS --fortran ${LEPTONFF_FORTRAN_EXAMPLE_COMPILER})
    endif()
    add_test(NAME installed_package COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_install.py ${CMAKE_CURRENT_BINARY_DIR}
      --cxx ${CMAKE_CXX_COMPILER} ${LEPTONFF_INSTALL_TEST_ARGUMENTS})
    add_test(NAME photonic_options COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_photonic_options.py)
    add_test(NAME independent_photonic_nll_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_photonic_nll_reference.py
      $<TARGET_FILE:leptonff_photonic_nll>)
    add_test(NAME independent_qcd_hadronic_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_qcd_hadronic_reference.py $<TARGET_FILE:leptonff_hadronic>)
    add_test(NAME independent_flavour_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_flavour_reference.py $<TARGET_FILE:leptonff_flavours>)
    add_test(NAME spectrum_analysis COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_spectrum_analysis.py)
    add_test(NAME comparison_analysis COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_comparison_analysis.py)
    add_test(NAME independent_lepton_spectrum_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_lepton_spectrum_reference.py
      $<TARGET_FILE:leptonff_lepton_spectrum> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf
      ${CMAKE_CURRENT_SOURCE_DIR}/configs/lepton_bins.conf
      --output ${CMAKE_CURRENT_BINARY_DIR}/lepton_spectrum_endpoint_checked_coarse.json)
    add_test(NAME independent_mass_spectrum_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_mass_spectrum_reference.py
      $<TARGET_FILE:leptonff_mass_spectrum> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf
      ${CMAKE_CURRENT_SOURCE_DIR}/configs/mass_bins.conf)
    add_test(NAME independent_coupled_decay_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_coupled_decay_reference.py
      $<TARGET_FILE:leptonff_coupled> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_pair_moment_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_pair_moment_reference.py
      $<TARGET_FILE:leptonff_pair_moments> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_photon_pair_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_photon_pair_reference.py
      $<TARGET_FILE:leptonff_photon_pair> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_xspace_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_xspace_reference.py
      $<TARGET_FILE:leptonff_xspace> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf --coarse-only)
    add_test(NAME independent_moment_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_moment_reference.py
      $<TARGET_FILE:leptonff_moments>)
    add_test(NAME independent_nlo_kernel_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_nlo_kernel_reference.py
      $<TARGET_FILE:leptonff_kernels>)
    add_test(NAME independent_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_reference.py
      $<TARGET_FILE:leptonff_decay> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_ll_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_ll_reference.py
      $<TARGET_FILE:leptonff_ll> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_matching_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_matching_reference.py
      $<TARGET_FILE:leptonff_match> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_vector_mass_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_vector_mass_reference.py
      $<TARGET_FILE:leptonff_hard> ${CMAKE_CURRENT_SOURCE_DIR}/configs/vector_mass.conf)
    add_test(NAME independent_fiducial_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_fiducial_reference.py
      $<TARGET_FILE:leptonff_fiducial> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME independent_hard_dressed_reference COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/python/check_hard_dressed_reference.py
      $<TARGET_FILE:leptonff_dressed> ${CMAKE_CURRENT_SOURCE_DIR}/configs/reference.conf
      --coarse-only)
  endif()
endif()
