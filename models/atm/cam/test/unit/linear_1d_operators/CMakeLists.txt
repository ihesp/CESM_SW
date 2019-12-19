# Local pFUnit files.
set(pf_sources
  test_diagonal.pf test_derivatives.pf test_arithmetic.pf)

# Sources to test.
set(sources_needed
  coords_1d.F90 linear_1d_operators.F90)
extract_sources("${sources_needed}" "${cam_sources}" test_sources)

# Do source preprocessing and add the executable.
create_pFUnit_test(linear_1d_operators linear_1d_operators_exe
  "${pf_sources}" "${test_sources}")

target_link_libraries(linear_1d_operators_exe csm_share)
