add_library(openggcm emfields.cxx)

target_include_directories(
  openggcm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
                  $<INSTALL_INTERFACE:include>)
