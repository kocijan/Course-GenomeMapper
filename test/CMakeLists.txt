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

include(GoogleTest)
gtest_discover_tests(empty_test)