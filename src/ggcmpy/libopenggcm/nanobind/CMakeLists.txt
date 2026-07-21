execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_openggcm bind.cxx)
target_link_libraries(_openggcm PRIVATE openggcm)

nanobind_add_stub(
  _openggcm_stub
  INSTALL_TIME
  MODULE
  _openggcm
  OUTPUT
  ggcmpy/_openggcm.pyi
  PYTHON_PATH
  $<TARGET_FILE_DIR:_openggcm>
  MARKER_FILE
  py.typed
  DEPENDS
  _openggcm)

install(TARGETS _openggcm DESTINATION ${SKBUILD_PROJECT_NAME})
