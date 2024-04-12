cmake_minimum_required(VERSION 3.15...3.27)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# F2PY headers
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
    OUTPUT _jrrlemodule.c _jrrle-f2pywrappers.f
    DEPENDS
        ggcmpy/backends/jrrle/_fortfile.F90
        ggcmpy/backends/jrrle/_jrrle.F90
    VERBATIM
    COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
            "${CMAKE_CURRENT_SOURCE_DIR}/ggcmpy/backends/jrrle/_fortfile.F90"
            "${CMAKE_CURRENT_SOURCE_DIR}/ggcmpy/backends/jrrle/_jrrle.F90"
            -m _jrrle --lower)

python_add_library(_jrrle MODULE
    "${CMAKE_CURRENT_BINARY_DIR}/_jrrlemodule.c"
    ggcmpy/backends/jrrle/_fortfile.F90
    ggcmpy/backends/jrrle/_jrrle.F90
    WITH_SOABI)
target_link_libraries(_jrrle PRIVATE fortranobject)

install(TARGETS _jrrle DESTINATION ${SKBUILD_PROJECT_NAME})
