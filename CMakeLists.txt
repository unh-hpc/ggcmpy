cmake_minimum_required(VERSION 3.15...3.27)
project(${SKBUILD_PROJECT_NAME} LANGUAGES Fortran)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

Python_add_library(_libggcm MODULE
    ggcmpy/test.F90
    ggcmpy/fortfile.F90
    ggcmpy/jrrle.F90
    WITH_SOABI)

install(TARGETS _libggcm DESTINATION ${SKBUILD_PROJECT_NAME})
