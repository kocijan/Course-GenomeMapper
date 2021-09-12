include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        e2239ee6043f73722e7aa812a459f54a28552929 # release-1.11.0
)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  empty_test
  empty_test.cpp
)
target_link_libraries(
  empty_test
  gtest_main
)

add_executable(
  alignment_test
  alignment_test.cpp
  ${PROJECT_SOURCE_DIR}/include/crimson_alignment_engine.hpp
)
target_link_libraries(
  alignment_test
  PUBLIC
  gtest_main
)

add_executable(
  minimizer_test
  minimizer_test.cpp
  ${PROJECT_SOURCE_DIR}/include/crimson_minimizer_engine.hpp
)
target_link_libraries(
  minimizer_test
  PUBLIC
  gtest_main
)

target_link_libraries(alignment_test PUBLIC crimson_alignment_engine)
target_link_libraries(minimizer_test PUBLIC crimson_minimizer_engine)

target_include_directories(alignment_test PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(minimizer_test PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_compile_options(alignment_test PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic -Wconversion>
)
target_compile_options(minimizer_test PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
  $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -pedantic -Wconversion>
)

include(GoogleTest)
gtest_discover_tests(empty_test)
gtest_discover_tests(alignment_test)
gtest_discover_tests(minimizer_test)
