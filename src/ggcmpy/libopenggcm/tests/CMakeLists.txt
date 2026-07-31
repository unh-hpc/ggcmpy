add_executable(test_emfields test_emfields.cxx)
target_link_libraries(test_emfields PRIVATE openggcm GTest::gtest_main)
gtest_discover_tests(test_emfields)

add_executable(test_tracing test_tracing.cxx)
target_link_libraries(test_tracing PRIVATE openggcm GTest::gtest_main)
gtest_discover_tests(test_tracing)
