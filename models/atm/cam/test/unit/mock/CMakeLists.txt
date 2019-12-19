string(REGEX REPLACE "[^;]+shr_sys_mod.F90"
  "${CMAKE_CURRENT_SOURCE_DIR}/shr_sys_mod.nompi_abortthrows.F90"
  share_sources "${share_sources}")

sourcelist_to_parent(share_sources)
