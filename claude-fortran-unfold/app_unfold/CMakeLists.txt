set(projectdir ${PROJECT_SOURCE_DIR})

set(sources-fpp
  initunfold.F90
  radial.F90
  unfold.F90)

set(fypp_flags ${FYPP_BUILD_FLAGS})
list(APPEND fypp_flags -I${projectdir}/src/dftbp/include)

dftbp_preprocess("${FYPP}" "${fypp_flags}" "F90" "f90" "${sources-fpp}" sources-f90-preproc)

add_executable(unfold ${sources-f90-preproc})

if(WITH_MPI)
  target_link_libraries(unfold PRIVATE MPI::MPI_Fortran)
endif()

target_link_libraries(unfold PRIVATE dftbplus)

install(TARGETS unfold DESTINATION "${CMAKE_INSTALL_BINDIR}")
