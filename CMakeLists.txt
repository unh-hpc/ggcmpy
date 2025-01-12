cmake_minimum_required(VERSION 3.15...3.27)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Ensure scikit-build modules
if(NOT SKBUILD)
  find_package(PythonInterp 3.10 REQUIRED)
  # Kanged -->
  # https://github.com/Kitware/torch_liberator/blob/master/CMakeLists.txt If
  # skbuild is not the driver; include its utilities in CMAKE_MODULE_PATH
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
    OUTPUT_VARIABLE SKBLD_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  list(APPEND CMAKE_MODULE_PATH "${SKBLD_DIR}/resources/cmake")
  message(STATUS "Looking in ${SKBLD_DIR}/resources/cmake for CMake modules")
endif()

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

# scikit-build style includes find_package(PythonExtensions REQUIRED) # for
# ${PYTHON_EXTENSION_MODULE_SUFFIX}

# Grab the variables from a local Python installation NumPy headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE)
# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_custom_command(
  OUTPUT _jrrlemodule.c _jrrle-f2pywrappers.f
  DEPENDS src/ggcmpy/backends/jrrle/_fortfile.F90
          src/ggcmpy/backends/jrrle/_jrrle.F90
  VERBATIM
  COMMAND
    "${Python_EXECUTABLE}" -m numpy.f2py
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ggcmpy/backends/jrrle/_fortfile.F90"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ggcmpy/backends/jrrle/_jrrle.F90" -m
    _jrrle --lower)

python_add_library(
  _jrrle
  MODULE
  ${CMAKE_CURRENT_BINARY_DIR}/_jrrlemodule.c
  src/ggcmpy/backends/jrrle/_fortfile.F90
  src/ggcmpy/backends/jrrle/_jrrle.F90
  ${F2PY_INCLUDE_DIR}/fortranobject.c
  WITH_SOABI)
target_include_directories(
  _jrrle PUBLIC ${F2PY_INCLUDE_DIR} ${NumPy_INCLUDE_DIRS}
                ${PYTHON_INCLUDE_DIRS})

install(TARGETS _jrrle DESTINATION ${SKBUILD_PROJECT_NAME})
