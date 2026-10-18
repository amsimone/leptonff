if(NOT LEPTONFF_BUILD_MATCHED)
  foreach(name IN ITEMS distributions fragmentation nlo_lepton nlo_qed moment_evolution flavours)
    add_executable(test_${name} test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE leptonff::core)
    add_test(NAME ${name} COMMAND test_${name})
  endforeach()
else()
  enable_language(C)
  add_executable(test_inclusive_mass_response test_inclusive_mass_response.cpp)
  target_link_libraries(test_inclusive_mass_response PRIVATE leptonff)
  add_test(NAME inclusive_mass_response COMMAND test_inclusive_mass_response)
  add_executable(test_c_api test_c_api.c)
  target_link_libraries(test_c_api PRIVATE leptonff)
  if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(test_c_api PRIVATE -UNDEBUG)
  endif()
  add_test(NAME c_api COMMAND test_c_api)
  foreach(name IN ITEMS response_table drell_yan production_records production_decay photonic_nll qcd_hadronic flavours flavour_evolution lab_spectrum fraction_cells coupled_lab photonic_dressing spin_decay lab_fragmentation endpoint_geometry lepton_spectrum mass_spectrum coupled_decay pair_moment photon_pair distributions fragmentation nlo_lepton nlo_qed moment_evolution xspace_evolution decay evolution matching vector_decay fiducial_vector hard_dressing)
    add_executable(test_${name} test_${name}.cpp)
    target_link_libraries(test_${name} PRIVATE leptonff)
    add_test(NAME ${name} COMMAND test_${name})
  endforeach()
  if(LEPTONFF_BUILD_REFERENCE_PROGRAMS)
    add_test(NAME spin_decay_reference COMMAND leptonff_spin)
    add_test(NAME nlo_kernel_reference COMMAND leptonff_kernels)
    add_test(NAME photon_pair_reference COMMAND leptonff_photon_pair
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME pair_moment_reference COMMAND leptonff_pair_moments
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME coupled_decay_reference COMMAND leptonff_coupled
      ${PROJECT_SOURCE_DIR}/configs/reference.conf 256 2 --dump)
    add_test(NAME lepton_spectrum_reference COMMAND leptonff_lepton_spectrum
      ${PROJECT_SOURCE_DIR}/configs/reference.conf ${PROJECT_SOURCE_DIR}/configs/lepton_bins.conf 256 1 --dump)
    add_test(NAME mass_spectrum_reference COMMAND leptonff_mass_spectrum
      ${PROJECT_SOURCE_DIR}/configs/reference.conf ${PROJECT_SOURCE_DIR}/configs/mass_bins.conf 256 1 --dump)
    add_test(NAME moment_reference COMMAND leptonff_moments)
    add_test(NAME xspace_reference COMMAND leptonff_xspace
      ${PROJECT_SOURCE_DIR}/configs/reference.conf 64)
    add_test(NAME decay_reference COMMAND leptonff_decay
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME ll_reference COMMAND leptonff_ll
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME matching_reference COMMAND leptonff_match
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME vector_mass_reference COMMAND leptonff_hard
      ${PROJECT_SOURCE_DIR}/configs/vector_mass.conf)
    add_test(NAME vector_mass_reject_angular_cuts COMMAND leptonff_hard
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    set_tests_properties(vector_mass_reject_angular_cuts PROPERTIES WILL_FAIL TRUE)
    add_test(NAME fiducial_reference COMMAND leptonff_fiducial
      ${PROJECT_SOURCE_DIR}/configs/reference.conf)
    add_test(NAME hard_dressed_reference COMMAND leptonff_dressed
      ${PROJECT_SOURCE_DIR}/configs/reference.conf 128 0.15 32 4)
  endif()
  find_package(Python3 3.10 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    if(LEPTONFF_BUILD_TOOLS)
      add_test(NAME production_replay_cli COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/tests/test_replay_cli.py $<TARGET_FILE:leptonff_replay>)
    endif()
    add_test(NAME dyturbo_source_preparation COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_dyturbo_prepare.py)
    add_test(NAME dyturbo_interface_contract COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_dyturbo_interface.py)
    add_test(NAME dyturbo_precision_analysis COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_dyturbo_precision.py)
    add_test(NAME fsr_reference_analysis COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_fsr_reference.py)
    find_program(LEPTONFF_FORTRAN_EXAMPLE_COMPILER NAMES gfortran-16 gfortran)
    set(LEPTONFF_INSTALL_TEST_ARGUMENTS "")
    if(LEPTONFF_FORTRAN_EXAMPLE_COMPILER)
      list(APPEND LEPTONFF_INSTALL_TEST_ARGUMENTS --fortran ${LEPTONFF_FORTRAN_EXAMPLE_COMPILER})
    endif()
    add_test(NAME installed_package COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_install.py ${PROJECT_BINARY_DIR}
      --cxx ${CMAKE_CXX_COMPILER} ${LEPTONFF_INSTALL_TEST_ARGUMENTS}
      --tools ${LEPTONFF_BUILD_TOOLS})
    add_test(NAME photonic_options COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_photonic_options.py)
    if(LEPTONFF_BUILD_REFERENCE_PROGRAMS)
      execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
        RESULT_VARIABLE LEPTONFF_AUDIT_PYTHON_STATUS OUTPUT_QUIET ERROR_QUIET)
      if(LEPTONFF_AUDIT_PYTHON_STATUS EQUAL 0)
        add_test(NAME radiator_validation_analysis COMMAND ${Python3_EXECUTABLE}
          ${PROJECT_SOURCE_DIR}/tests/test_radiator_validation.py)
        add_test(NAME mcmule_fixed_order_analysis COMMAND ${Python3_EXECUTABLE}
          ${PROJECT_SOURCE_DIR}/tests/test_mcmule_analysis.py)
        add_test(NAME independent_complete_mass_radiator COMMAND ${Python3_EXECUTABLE}
          ${PROJECT_SOURCE_DIR}/adapters/validation/check_mass_radiator.py
          --client $<TARGET_FILE:leptonff_mass_audit> --grids 4096 8192)
      endif()
      add_test(NAME independent_photonic_nll_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_photonic_nll_reference.py
        $<TARGET_FILE:leptonff_photonic_nll>)
      add_test(NAME independent_qcd_hadronic_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_qcd_hadronic_reference.py $<TARGET_FILE:leptonff_hadronic>)
      add_test(NAME independent_flavour_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_flavour_reference.py $<TARGET_FILE:leptonff_flavours>)
    endif()
    add_test(NAME spectrum_analysis COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_spectrum_analysis.py)
    add_test(NAME comparison_analysis COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_comparison_analysis.py)
    if(LEPTONFF_BUILD_REFERENCE_PROGRAMS)
      add_test(NAME independent_lepton_spectrum_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_lepton_spectrum_reference.py
        $<TARGET_FILE:leptonff_lepton_spectrum> ${PROJECT_SOURCE_DIR}/configs/reference.conf
        ${PROJECT_SOURCE_DIR}/configs/lepton_bins.conf
        --output ${PROJECT_BINARY_DIR}/lepton_spectrum_endpoint_checked_coarse.json)
      add_test(NAME independent_mass_spectrum_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_mass_spectrum_reference.py
        $<TARGET_FILE:leptonff_mass_spectrum> ${PROJECT_SOURCE_DIR}/configs/reference.conf
        ${PROJECT_SOURCE_DIR}/configs/mass_bins.conf)
      add_test(NAME independent_coupled_decay_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_coupled_decay_reference.py
        $<TARGET_FILE:leptonff_coupled> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_pair_moment_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_pair_moment_reference.py
        $<TARGET_FILE:leptonff_pair_moments> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_photon_pair_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_photon_pair_reference.py
        $<TARGET_FILE:leptonff_photon_pair> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_xspace_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_xspace_reference.py
        $<TARGET_FILE:leptonff_xspace> ${PROJECT_SOURCE_DIR}/configs/reference.conf --coarse-only)
      add_test(NAME independent_moment_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_moment_reference.py
        $<TARGET_FILE:leptonff_moments>)
      add_test(NAME independent_nlo_kernel_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_nlo_kernel_reference.py
        $<TARGET_FILE:leptonff_kernels>)
      add_test(NAME independent_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_reference.py
        $<TARGET_FILE:leptonff_decay> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_ll_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_ll_reference.py
        $<TARGET_FILE:leptonff_ll> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_matching_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_matching_reference.py
        $<TARGET_FILE:leptonff_match> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_vector_mass_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_vector_mass_reference.py
        $<TARGET_FILE:leptonff_hard> ${PROJECT_SOURCE_DIR}/configs/vector_mass.conf)
      add_test(NAME independent_fiducial_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_fiducial_reference.py
        $<TARGET_FILE:leptonff_fiducial> ${PROJECT_SOURCE_DIR}/configs/reference.conf)
      add_test(NAME independent_hard_dressed_reference COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/python/check_hard_dressed_reference.py
        $<TARGET_FILE:leptonff_dressed> ${PROJECT_SOURCE_DIR}/configs/reference.conf
        --coarse-only)
    endif()
  endif()
endif()

foreach(name IN ITEMS theory core_api core_evolution)
  add_executable(test_${name} test_${name}.cpp)
  target_link_libraries(test_${name} PRIVATE leptonff::core)
  add_test(NAME ${name} COMMAND test_${name})
endforeach()
target_compile_definitions(test_core_api PRIVATE LEPTONFF_EXPECTED_VERSION="${PROJECT_VERSION}")
find_package(Python3 3.10 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  if(NOT LEPTONFF_BUILD_MATCHED)
    add_test(NAME installed_core_package COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_install.py ${PROJECT_BINARY_DIR}
      --cxx ${CMAKE_CXX_COMPILER} --core-only --tools OFF)
  endif()
  add_test(NAME canonical_benchmark COMMAND ${Python3_EXECUTABLE}
    ${PROJECT_SOURCE_DIR}/tests/check_benchmark.py $<TARGET_FILE:leptonff_benchmark>)
  if(BUILD_SHARED_LIBS)
    add_test(NAME python_bindings COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_bindings.py $<TARGET_FILE:leptonff_core>)
    if(LEPTONFF_BUILD_MATCHED)
      add_test(NAME python_decay_bindings COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/tests/test_bindings.py $<TARGET_FILE:leptonff> --decay)
    endif()
  endif()
endif()
